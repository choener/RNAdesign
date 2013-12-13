{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

module BioInf.RNAdesign where

import qualified Data.Array.IArray as A
import System.IO.Unsafe
import Control.Monad.IO.Class
import Control.Monad.Primitive.Class
import Control.Monad.Primitive
import System.Random.MWC.Monad
import Control.Monad
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Data.List (sort,group)
import qualified Data.Map as M
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import Control.Arrow
import System.IO.Unsafe -- TODO remove
import Data.List
import Data.Tuple.Select

import Biobase.Primary
import Biobase.Secondary.Diagrams
import Biobase.Secondary
import Biobase.Vienna
import qualified BioInf.ViennaRNA.Bindings  as RNA     -- NOTE removes the ability to call into ghci!
import BioInf.ViennaRNA.Eval

import BioInf.RNAdesign.Graph
import BioInf.RNAdesign.OptParser
import BioInf.RNAdesign.Assignment
import BioInf.RNAdesign.LogMultinomial

import Debug.Trace



-- A single candidate, with its sequence and the score, this sequence receives.
-- Candidates are ordered by their scores.

data Candidate = Candidate
  { candidate :: Primary
  , score :: Score
  } deriving (Eq,Show)

instance Ord Candidate where
  (Candidate _ a) <= (Candidate _ b) = ropt a <= ropt b

-- | Create an initial, legal, candidate. Give it a really bad score.

mkInitial :: (MonadPrim m, PrimMonad m) => Int -> DesignProblem -> Rand m Candidate
mkInitial l dp = do
  let z = VU.replicate l nA
  c <- foldM mutateOneAssignment z $ assignments dp
  return $ Candidate c (Score [] 999999)

{-
-- | Sum probabilities over base pairs in the structural constraints

sumProbStructures :: Primary -> [D1Secondary] -> Double
sumProbStructures inp ss = s where
  s = sum $ map ((bp A.!) . first (+1) . second (+1)) ps
  ps = concatMap snd (map fromD1S ss :: [(Int,[PairIdx])])
  bp = let (_,_,bp') = unsafePerformIO (RNA.part $ concatMap show $ VU.toList inp) in bp'

sumProbNotStructures :: Primary -> [D1Secondary] -> Double
sumProbNotStructures inp ss = undefined

probabilityDefect inp str = s where
  s = sum (map (bp A.!) ps) - sum (map (bp A.!) ups)
  ups = [ (i,j) | i<-[1..l], j<-[i..l] ] \\ ps
  (l,ps) = second (map (first (+1) . second (+1))) $ fromD1S str :: (Int,[PairIdx])
  bp = let (_,_,bp') = unsafePerformIO (RNA.part $ concatMap show $ VU.toList inp) in bp'
-}

probabilityDefectAll inp ss = s where
  ca :: A.Array (Int,Int) Double
  ca = A.amap (\c -> c / n) . A.accumArray (+) 0 ((1,1),(l,l)) $ zip ps (repeat 1)
  n = genericLength ss
--  s = sum (map (abs . (n-) . (bp A.!)) ps) + sum (map (bp A.!) ups)
  s = sum (map (\ix -> abs $ ca A.! ix - bp A.! ix) ps) + sum (map (bp A.!) ups)
  l = VU.length inp
  ups = [ (i,j) | i<-[1..l], j<-[i..l] ] \\ ps
  ps = map (first (+1) . second (+1)) $ concatMap snd (map fromD1S ss :: [(Int,[PairIdx])])
  bp = let (_,_,bp') = unsafePerformIO (RNA.part $ concatMap show $ VU.toList inp) in bp'

ensembleDefect inp str = s where
  s = n - 2 * sps - sus
  n = fromIntegral $ VU.length inp
  sps = sum $ map (bp A.!) ps
  sus = sum $ [bp A.! (i,j) | i <- us, j <- [i..n]]
  ps = map (first (+1) . second (+1)) $ snd $ (fromD1S str :: (Int,[PairIdx]))
  us = [1..n] \\ (map fst ps ++ map snd ps)
  (_,_,bp) = unsafePerformIO (RNA.part $ concatMap show $ VU.toList inp)

-- | Resolve the optimization task. Each possible optimization function is
-- given here. Try to keep the functions defined here in sync with some
-- (non-existent ;-) documentation.

resolveOpt :: String -> t -> Primary -> [D1Secondary] -> Double
resolveOpt optfun ener inp secs = parseOptString l sops mops gops props optfun where
  l = length secs
  sops =
    [ ("eos"   , \k -> unsafePerformIO $ RNA.eos (concatMap show (VU.toList inp)) (fromD1S $ secs !! (k-1)))
    , ("ed"    , \k -> ensembleDefect inp (secs !! (k-1))) -- ensemble defect
--    , ("pdef"  , \k -> probabilityDefect inp (secs !! (k-1)))
--    [ ("EOS",\k -> let (Deka e) = fst $ rnaEval ener inp $ secs !! (k-1) in fromIntegral e / 100)
--    , ("PF" ...
    ]
  mops =
    [ ("sum",sum)
    , ("max",maximum)
    , ("min",minimum)
    ]
  gops =
    [ ("Ged"   , probabilityDefectAll inp secs) -- global ensemble defect a la ``me''
    , ("gibbs" , sel1 . unsafePerformIO $ RNA.part (concatMap show (VU.toList inp)))
    , ("mfe"   , fst  . unsafePerformIO $ RNA.mfe (concatMap show (VU.toList inp)))
--    , ("Pin" , sumProbStructures inp secs)
    ]
  props =
    [ ("logMN", \ps -> lmn ps inp)
    ]

lmn ps inp = logMultinomial l p c where
  l   = VU.length inp
  p   = VU.fromList ps
  cM  = M.fromList . map (\z -> (head z, length z)) . group . sort $ VU.toList inp
  c   = VU.fromList $ map (\z -> M.findWithDefault 0 z cM) acgu

data Score = Score
  { eoss :: [Deka]
  , ropt :: Double
  } deriving (Eq,Show,Read)

instance Ord Score where
  (Score _ a) <= (Score _ b) = a<=b

scoreSequence :: String -> Vienna2004 -> DesignProblem -> Primary -> Score
scoreSequence optfun ener DesignProblem{..} s = score where
  score = Score
    { eoss = error "don't call this" -- map (fst . rnaEval ener s) structures
    , ropt = resolveOpt optfun ener s structures
    }

-- | This structure defines a "design problem"

data DesignProblem = DesignProblem
  { structures  :: [D1Secondary]
  , assignments :: [Assignment]
  } deriving (Eq,Read,Show)

-- | Given a set of structures, create the set of independent graphs and
-- assignment possibilities.

mkDesignProblem :: Int -> [String] -> [String] -> DesignProblem
mkDesignProblem asnLimit xs scs = dp where
  dp = DesignProblem
        { structures  = map mkD1S xs
        , assignments = as
        }
  gs = independentGraphs xs
  as = map (allCandidates asnLimit sv) gs
  ss = M.map fixup . M.unionsWith ((nub .) . (++)) $ map (M.fromList . zip [0..] . (map ((:[]). mkNuc))) scs
  sv = V.fromList $ map (\k -> M.findWithDefault acgu k ss) [0 .. length (head xs) - 1]
  fixup zs = filter (/=nN) $ if (all (==nN) zs) then acgu else zs

unfoldStreamNew
  :: forall m . (MonadPrim m, PrimMonad m)
  => Int -> Int -> Int -> (Primary -> Score) -> (Candidate -> Candidate -> Rand m Bool) -> DesignProblem -> Candidate -> SM.Stream (Rand m) Candidate
unfoldStreamNew burnin number thin score f dp = go where
  go s  = SM.map snd                                                          -- remove remaining indices from stream
        . SM.take number                                                      -- take the number of sequences we want
        . SM.filter ((==0) . flip mod thin . fst)                             -- keep only every thin'th sequence
        . SM.indexed                                                          -- add index
        . SM.drop burnin                                                      -- drop the burnin sequences
        . SM.drop 1                                                           -- drop original input
        . SM.scanlM' (mutateOneAssignmentCandidateWith score f) s             -- starting with 's', mutate s further and further using cycled assignments
        $ SM.unfoldr (Just . first head . splitAt 1) (cycle $ assignments dp) -- create inifinite cycled assignments

unfoldStream
  :: forall m . (MonadPrim m, PrimMonad m)
  => Int -> Int -> Int -> (Primary -> Primary -> Rand m Bool) -> DesignProblem -> Primary -> SM.Stream (Rand m) Primary
unfoldStream burnin number thin f dp = go where
  go s  = SM.map snd                                                          -- remove remaining indices from stream
        . SM.take number                                                      -- take the number of sequences we want
        . SM.filter ((==0) . flip mod thin . fst)                             -- keep only every thin'th sequence
        . SM.indexed                                                          -- add index
        . SM.drop burnin                                                      -- drop the burnin sequences
        . SM.drop 1                                                           -- drop original input
        . SM.scanlM' (mutateOneAssignmentWith f) s                            -- starting with 's', mutate s further and further using cycled assignments
        $ SM.unfoldr (Just . first head . splitAt 1) (cycle $ assignments dp) -- create inifinite cycled assignments

-- | Mutate the sequence in a candidate

mutateOneAssignmentCandidateWith
  :: (MonadPrim m, PrimMonad m)
  => (Primary -> Score) -> (Candidate -> Candidate -> Rand m Bool) -> Candidate -> Assignment -> Rand m Candidate
mutateOneAssignmentCandidateWith score f old Assignment{..} = do
  i <- uniformR (0,V.length assignment -1) -- inclusive range for Int
  let cs = VU.zip columns (assignment V.! i)
  let nw = VU.update (candidate old) cs
  let new = Candidate nw (score nw)
  b <- f old new
  return $ if b then new else old

-- | Mutate the sequence using one assignment with evaluation function.

mutateOneAssignmentWith
  :: (MonadPrim m, PrimMonad m)
  => (Primary -> Primary -> Rand m Bool) -> Primary -> Assignment -> Rand m Primary
mutateOneAssignmentWith f old Assignment{..} = do
  i <- uniformR (0,V.length assignment -1) -- inclusive range for Int
  let cs = VU.zip columns (assignment V.! i)
  let new = VU.update old cs
  b <- f old new
  return $ if b then new else old

-- | Create a number of sequences, thinning the list of candidates to yield
-- more independent candidates. The optimization function is used to make the
-- choice between emitting the current candidate again and selecting a new one.

generateSequences
  :: (MonadPrim m, PrimMonad m)
  => Int -> Int -> (Primary -> Primary -> Rand m Bool) -> DesignProblem -> Primary -> Rand m [Primary]
generateSequences number thin f dp s = go number thin s where
  go n t s
    | n < 1 = return []
    | t == 0 = do s' <- mutateSequence f dp s
                  ss <- go (n-1) thin s'
                  return $ s' : ss
    | otherwise = mutateSequence f dp s >>= go n (t-1)

-- | Mutate a sequence using the possible assignments.

mutateSequence
  :: (MonadPrim m, PrimMonad m)
  => (Primary -> Primary -> Rand m Bool) -> DesignProblem -> Primary -> Rand m Primary
mutateSequence f dp old = do
  new <- foldM mutateOneAssignment old $ assignments dp
  b <- f old new
  return $ if b then new else old

-- | Mutate the sequence using one assignment.

mutateOneAssignment
  :: (MonadPrim m, PrimMonad m)
  => Primary -> Assignment -> Rand m Primary
mutateOneAssignment s Assignment{..} = do
  i <- uniformR (0,V.length assignment -1) -- inclusive range for Int
  let cs = VU.zip columns (assignment V.! i)
  return $ VU.update s cs


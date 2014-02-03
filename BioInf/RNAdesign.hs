{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

module BioInf.RNAdesign where

import           Control.Arrow (first,second)
import           Control.Monad.Primitive
import           Control.Monad.Primitive.Class
import           Data.List (nub,group,sort,(\\),genericLength)
import           Data.Tuple.Select (sel1)
import qualified Data.Array.IArray as A
import qualified Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import           System.IO.Unsafe (unsafePerformIO)
import           System.Random.MWC.Monad

import           Biobase.Primary
import           Biobase.Primary.IUPAC
import           Biobase.Secondary.Diagrams
import           Biobase.Secondary (PairIdx(..))
import           Biobase.Vienna
import qualified BioInf.ViennaRNA.Bindings  as RNA

import           BioInf.RNAdesign.Assignment
import           BioInf.RNAdesign.CandidateChain
import           BioInf.RNAdesign.Graph
import           BioInf.RNAdesign.LogMultinomial
import           BioInf.RNAdesign.OptParser



-- |

probabilityDefectAll inp ss = s where
  ca :: A.Array (Int,Int) Double
  ca = A.amap (\c -> c / n) . A.accumArray (+) 0 ((1,1),(l,l)) $ zip ps (repeat 1)
  n = genericLength ss
  s = sum (map (\ix -> abs $ ca A.! ix - bp A.! ix) ps) + sum (map (bp A.!) ups)
  l = VU.length inp
  ups = [ (i,j) | i<-[1..l], j<-[i..l] ] \\ ps
  ps = map (first (+1) . second (+1)) $ concatMap snd (map fromD1S ss :: [(Int,[PairIdx])])
  bp = let (_,_,bp') = unsafePerformIO (RNA.part $ concatMap show $ VU.toList inp) in bp'

-- |

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
    ]
  props =
    [ ("logMN", \ps -> lmn ps inp)
    ]

-- |

lmn ps inp = logMultinomial l p c where
  l   = VU.length inp
  p   = VU.fromList ps
  cM  = M.fromList . map (\z -> (head z, length z)) . group . sort $ VU.toList inp
  c   = VU.fromList $ map (\z -> M.findWithDefault 0 z cM) acgu

-- |

scoreSequence :: String -> Vienna2004 -> DesignProblem -> Primary -> Score
scoreSequence optfun ener DesignProblem{..} s = score where
  score = Score $ resolveOpt optfun ener s structures

-- | Given a set of structures, create the set of independent graphs and
-- assignment possibilities.

mkDesignProblem :: Int -> [String] -> String -> DesignProblem
mkDesignProblem asnLimit xs scs = dp where
  dp = DesignProblem
        { structures  = map mkD1S xs
        , assignments = as
        }
  gs = independentGraphs xs
  as = map (allCandidates asnLimit sv) gs
  --ss = M.map fixup . M.unionsWith ((nub .) . (++)) $ map (M.fromList . zip [0..] . (map ((:[]). mkNuc))) scs
  ss = M.map fixup . M.fromList . zip [0..] . map (map mkNuc . fromSymbol) $ scs
  sv = V.fromList $ map (\k -> M.findWithDefault acgu k ss) [0 .. length (head xs) - 1]
  fixup zs = filter (/=nN) $ if (all (==nN) zs) then acgu else zs


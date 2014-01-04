{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RecordWildCards #-}

module BioInf.RNAdesign.CandidateChain where

import           Control.Arrow (first)
import           Control.Monad (foldM)
import           Control.Monad.Primitive
import           Control.Monad.Primitive.Class
import           Data.Function (on)
import qualified Data.Vector as V
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import qualified Data.Vector.Unboxed as VU
import           System.Random.MWC.Monad

import           Biobase.Primary
import           Biobase.Secondary.Diagrams
import           Biobase.Vienna

import           BioInf.RNAdesign.Assignment (Assignment(..))



-- | A single candidate, with its sequence and the score, this sequence
-- receives.  Candidates are ordered by their scores.

data Candidate = Candidate
  { candidate :: Primary
  , score :: Score
  } deriving (Eq,Show)

instance Ord Candidate where
  (<=) = (<=) `on` score

-- | The likelihood score we get.
--
-- TODO replace Score Likelihood / LogLikelihood (once we switch to the more
-- generic MCMC library)

newtype Score = Score { unScore :: Double }
  deriving (Eq,Ord,Show,Read)

-- | This structure defines a "design problem"

data DesignProblem = DesignProblem
  { structures  :: [D1Secondary]
  , assignments :: [Assignment]
  } deriving (Eq,Read,Show)

-- | Create an initial, legal, candidate.

mkInitial :: (MonadPrim m, PrimMonad m) => (Primary -> Score) -> Int -> DesignProblem -> Rand m Candidate
mkInitial scoring l dp = do
  let z = VU.replicate l nA
  foldM (mutateOneAssignmentCandidateWith scoring (\_ _ -> return True)) (Candidate z (scoring z)) $ assignments dp

-- | Create a stream of 'Candidate's from an initial candidate.

unfoldStreamCandidate
  :: forall m . (MonadPrim m, PrimMonad m)
  => Int -> Int -> Int -> (Primary -> Score) -> (Candidate -> Candidate -> Rand m Bool) -> DesignProblem -> Candidate
  -> SM.Stream (Rand m) Candidate
unfoldStreamCandidate burnin number thin score f dp = go where
  go s  = SM.map snd                                                          -- remove remaining indices from stream
        . SM.take number                                                      -- take the number of sequences we want
        . SM.filter ((==0) . flip mod thin . fst)                             -- keep only every thin'th sequence
        . SM.indexed                                                          -- add index
        . SM.drop burnin                                                      -- drop the burnin sequences
        . SM.drop 1                                                           -- drop original input
        . SM.scanlM' (mutateOneAssignmentCandidateWith score f) s             -- starting with 's', mutate s further and further using cycled assignments
        $ SM.unfoldr (Just . first head . splitAt 1) (cycle $ assignments dp) -- create inifinite cycled assignments

-- | Mutate the current (or "old") sequence under the possible 'assignment's as
-- prescribed by 'Assignment'. The modifying assignment is selected uniformly.
-- The monadic @old -> new -> Rand m Bool@ function chooses between the old and
-- the new candidate. It can be used to, e.g., allow always choosing "new" if
-- it is better, but choosing "new" as well if some stochastic value (hence
-- dependence on @Rand m@) indicates so.

mutateOneAssignmentCandidateWith
  :: (MonadPrim m, PrimMonad m)
  => (Primary -> Score)                       -- ^ the likelihood function, gives a sequence a score
  -> (Candidate -> Candidate -> Rand m Bool)  -- ^ choose between old and new sequence (monadic, stochastic)
  -> Candidate                                -- ^ "old" / current sequence
  -> Assignment                               -- ^ possible assignments for the sequence
  -> Rand m Candidate                         -- ^ the "new" sequence
mutateOneAssignmentCandidateWith score f old Assignment{..} = do
  i <- uniformR (0,V.length assignment -1) -- inclusive range for Int
  let cs = VU.zip columns (assignment V.! i)
  let nw = VU.update (candidate old) cs
  let new = Candidate nw (score nw)
  b <- f old new
  return $ if b then new else old


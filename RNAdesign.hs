{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

-- | Given one or more structural constraints, and possibly sequence
-- constraints for certain columns, design a sequence which is optimal
-- according to a user-defined optimization function. Optimization works via a
-- Markov Chain.

module Main where

import           Control.Arrow
import           Control.Monad
import           Data.Char (isAlpha)
import           Data.Function
import           Data.List
import           Data.Ord
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import qualified Data.Vector.Unboxed as VU
import           System.Console.CmdArgs
import           System.IO
import           System.Random.MWC.Distributions.Monad
import           System.Random.MWC.Monad
import           Text.Printf

import           Biobase.Primary
import           Biobase.Vienna
import qualified Biobase.Turner.Import as TI

import BioInf.RNAdesign
import BioInf.RNAdesign.CandidateChain
import BioInf.RNAdesign.Assignment



-- * Configuration

data Config = Config
  { number      :: Int
  , thin        :: Int
  , burnin      :: Int
  , scale       :: Double
  , optfun      :: String
  , veclen      :: Int
  , turner      :: String
--  , exhaustive  :: Bool     -- TODO want to think about this for number of structures > 3, IF the total sequence space size is less than say 100.000
  , initial     :: String
  , sequenceConstraints :: Bool
  , explore     :: Bool
  } deriving (Show,Data,Typeable)

config = Config
  { number      =  50 &= help "Number of candidate sequences to generate (50)"
  , thin        =  50 &= help "keep only every n'th sequence (50)"
  , burnin      = 100 &= help "remove the first burnin sequences (100)"
  , scale       =   1 &= help "acceptance scale parameter (1); exponentially distributed with mean 'scale^(-1)' (smaller scale means longer jumps)"
  , optfun      = "sum(eos,all)" &= help "optimization function, \"sum(eos,all)\" tries to minimize the sum of the energies"
  , veclen      = 1000000 &= help "multiple structure constraints lead to large connected components, veclen restricts the number of component solutions to store."
  , turner      = "./params" &= help "directory containing the Turner 2004 RNA energy tables (with a default of \"./params/\""
--  , exhaustive  = False &= help "exhaustively search the nucleotide space"
  , initial     = "" &= help "start from this initial sequence"
  , explore     = def &= help "explore sequences, do not sort of nub list"
  , sequenceConstraints = def &= help "activate sequence constraints"
  } &= help shortHelp &= details longHelp &= summary "RNAdesign v0.0.2, (C) Christian Hoener zu Siederdissen 2013, choener@tb.univie.ac.at" &= program "RNAdesign"

shortHelp = "The defaults work acceptably well and produce a results extremely fast. "

longHelp =
  [ "RNAdesign designs RNA sequences given one or more structural targets. The"
  , "program offers a variety of optimization functions that each can be used to"
  , "optimize candidate sequence towards a certain goal, say, minimal ensemble"
  , "defect or small energetic distance to another target structure. By giving a"
  , "complex \"--optun\", many different design goals can be tried. The following"
  , "functions are available:"
  , "binary, combining:"
  , "+ - * /  :: the four basic operations"
  , "^        :: (^) generalized power function"
  , ""
  , "binary, apply function to many targets:"
  , "sum max min   :: run function over set of targets: sum(eos,1,2) or sum(eos,all)"
  , ""
  , "unary, apply to single target:"
  , "eos      :: energy of a structure: eos(1)"
  , "ed       :: ensemble defect of a structure: ed(3)"
  , "nullary, constant for the current sequence:"
  , "Ged      :: global, weighted ensemble defect: Ged"
  , "gibbs    :: gibbs free energy of sequence"
  , "mfe      :: minimum free energy of sequence"
  , ""
  , "special:"
  , "logMN    :: requires four parameters logMN(0.2,0.3,0.3,0.2) penalizes according to given mono-nucleotide distribution"
  , ""
  , "A good optimization goal is (as an example for three targets):"
  , "--optfun \"eos(1)+eos(2)+eos(3) - 3*gibbs +"
  , "           1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)\""
  , "This way, the sequence produces close-to-mfe foldings with the targets (left) and the targets are close together in terms of energy. (1 *) scales the two terms according to user choice."
  , "\n\n\n"
  , "If you find this tool useful, please cite:"
  , ""
  , "Christian Hoener zu Siederdissen, Stefan Hammer, Ingrid Abfalter, Ivo L. Hofacker, Christoph Flamm, Peter F. Stadler."
  , "A Graph Coloring Approach to Designing Multi-Stable Nucleic Acid Sequences."
  , "submitted, 2013."
  , ""
  , "Contact: choener@tbi.univie.ac.at"
  , "Given one or more structures in dot-bracket format of the same length, returns a compatible assignment of nucleotides."
  , "Compatible nucleotides are those that allow folding of the sequence into all given structures."
  ]

main = do
  hSetBuffering stdout NoBuffering
  hSetBuffering stderr NoBuffering
  cmds@Config{..} <- cmdArgs config
  turner <- fmap turnerToVienna $ TI.fromDir turner "" ".dat" 
  strs' <- fmap lines $ getContents
  let (scs,strs) = partition (any isAlpha) . filter ((">"/=) . take 1) $ strs'
  unless (length strs > 0) $ error "no structures given!"
  let l = length $ head strs
  unless (all ((l==) . length) strs) $ error "structures of different size detected"
  let dp = mkDesignProblem veclen strs (if sequenceConstraints then scs else [])
  let defOpt old new = let oldS = scoreSequence optfun turner dp old
                           newS = scoreSequence optfun turner dp new
                       in  do t <- exponential scale
                              return $ unScore newS <= unScore oldS || t >= unScore newS - unScore oldS
  let calcScore = scoreSequence optfun turner dp
  let walk old new = do t <- exponential scale
                        let sn = unScore $ score new
                        let so = unScore $ score old
                        return $ sn <= so || t >= sn - so
  let ini = if null initial         -- start from initial sequence or generate one from the ensemble
              then mkInitial calcScore l dp
              else let pri = mkPrimary initial
                   in  return $ Candidate pri (calcScore pri)
  xs <- runWithSystemRandom . asRandIO $ (ini >>= SM.toList . unfoldStream burnin number thin calcScore walk dp)
  let pna = product . map numAssignments $ assignments dp
  printf "# Size of sequence space: %d %s\n\n" pna (show . map numAssignments $ assignments dp)
  unless (pna>0) $ error "empty sequence space, aborting!"
  mapM_ (\ys -> printf "%s %4d %8.2f\n" (concatMap show . VU.toList . candidate . head $ ys) (length ys) (unScore . score $ head ys))
    . ( if   explore
        then map (:[])
        else ( sortBy (comparing (unScore . score . head))
             . groupBy ((==) `on` candidate)
             . sortBy (comparing candidate)
             )
      )
    $ xs


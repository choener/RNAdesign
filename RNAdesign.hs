{-# LANGUAGE DoAndIfThenElse #-}
{-# LANGUAGE TemplateHaskell #-}
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
import           Data.FileEmbed
import           Data.Function
import           Data.List
import           Data.Ord
import           Data.Version (showVersion)
import qualified Data.ByteString.Char8 as BS
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

import Paths_RNAdesign (version)



-- * Configuration

data Config = Config
  { number              :: Int
  , thin                :: Int
  , burnin              :: Int
  , scale               :: Double
  , optfun              :: String
  , veclen              :: Int
  , turner              :: String
  , initial             :: String
  , sequenceConstraints :: Bool
  , explore             :: Bool
  , showManual          :: Bool
  } deriving (Show,Data,Typeable)

config = Config
  { number              =  50             &= help "Number of candidate sequences to generate (50)"
  , thin                =  50             &= help "keep only every n'th sequence (50)"
  , burnin              = 100             &= help "remove the first burnin sequences (100)"
  , scale               =   1             &= help "acceptance scale parameter (1); exponentially distributed with mean 'scale^(-1)' (smaller scale means longer jumps)"
  , optfun              = "sum(eos,all)"  &= help "optimization function, \"sum(eos,all)\" tries to minimize the sum of the energies"
  , veclen              = 1000000         &= help "multiple structure constraints lead to large connected components, veclen restricts the number of component solutions to store."
  , turner              = "./params"      &= help "directory containing the Turner 2004 RNA energy tables (with a default of \"./params/\""
  , initial             = ""              &= help "start from this initial sequence"
  , explore             = def             &= help "explore sequences, do not sort of nub list"
  , sequenceConstraints = def             &= help "activate sequence constraints"
  , showManual          = def             &= help ""
  } &= help shortHelp
    &= details [] -- longHelp
    &= summary ("RNAdesign " ++ showVersion version ++ " (C) Christian Hoener zu Siederdissen 2013--2014, choener@tbi.univie.ac.at")
    &= program "RNAdesign"

shortHelp = "The defaults work acceptably well and produce a results extremely fast. "

embeddedManual = $(embedFile "README.md")

main = do
  hSetBuffering stdout NoBuffering
  hSetBuffering stderr NoBuffering
  cmds@Config{..} <- cmdArgs config
  if showManual
  then BS.putStrLn embeddedManual
  else do
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


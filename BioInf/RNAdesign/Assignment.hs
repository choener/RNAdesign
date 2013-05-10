{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE ScopedTypeVariables #-}

module BioInf.RNAdesign.Assignment where

import Control.Arrow
import Data.Graph.Inductive.Graph
import Data.List (nub,sortBy,sort,genericLength)
import Data.Ord
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import Control.Lens.Tuple
import Control.Lens
import Data.Function

import Biobase.Primary
import Biobase.Secondary.Vienna
import Data.Graph.Inductive.Query.Ear

import BioInf.RNAdesign.Graph

import Data.Graph.Inductive.Query
import Debug.Trace



data Assignment = Assignment
  { columns        :: VU.Vector Int
  , assignment     :: V.Vector (VU.Vector Nuc)
  , isExhaustive   :: Bool
  , numAssignments :: Integer
  } deriving (Eq,Read,Show)

-- | Given a graph with base pairing constraints, return a 'Assignments' data
-- structure that provides all legal assignments.

allCandidates :: (DynGraph gr) => Int -> V.Vector [Nuc] -> gr () () -> Assignment
allCandidates maxC sv g
  | noNodes g == 1 = let [n] = nodes g
                         svn = sv V.! n
                     in  Assignment (VU.singleton n) (V.fromList $ map VU.singleton svn) True (genericLength svn)
  | noNodes g == 2 = let [n,m] = nodes g
                         svnm = [[a,b] | a <- sv V.! n, b <- sv V.! m, isViennaPair a b]
                     in  Assignment (VU.fromList [n,m]) (V.fromList $ map VU.fromList $ svnm) True (genericLength svnm)
  | noNodes g >  2 = let es :: [(Int,Int)] = nub . filter (uncurry (<)) . map (\(a,b,_) -> (a,b)) . sortBy (compare `on` sel3) . labEdges $ g' -- $ ears g
                         g' = case bcc g of
                                [_] -> ears g
                                xs  -> emap (const 0) g
                         o = sort $ mkEL es
                         (as,num) = second genericLength . splitAt maxC . map VU.fromList $ mkAssignments sv es
                     in  Assignment (VU.fromList o) (V.fromList $ take maxC as) (num==1) (genericLength as + num)

mkEL = nub . concatMap (\(a,b) -> [a,b])

mkAssignments sv es = map (map snd . sortBy (comparing fst)) $ mkA [] es where
  mkA :: [(Int,Nuc)] -> [(Int,Int)] -> [[(Int,Nuc)]]
  mkA dones [] = [[]]
  mkA dones ((a,b):ds)
    | Nothing <- a', Nothing <- b'
    = [ (a,n):(b,m):ns | n <- sv V.! a, m <- sv V.! b, isViennaPair n m, ns <- mkA ((a,n):(b,m):dones) ds ]
    | Nothing <- a', Just m <- b'
    = [ (a,n):ns | n <- sv V.! a, isViennaPair n m, ns <- mkA ((a,n):dones) ds ]
    | Just n <- a', Nothing <- b'
    = [ (b,m):ns | m <- sv V.! b, isViennaPair n m, ns <- mkA ((b,m):dones) ds ]
    | Just n <- a', Just m <- b'
    = if isViennaPair n m then mkA dones ds else []
    where
      a' = lookup a dones
      b' = lookup b dones

vps = filter (uncurry isViennaPair) [(a,b) | a<-cgau, b<-cgau]

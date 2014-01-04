
module BioInf.RNAdesign.Graph where

import Control.Arrow (first,second)
import Data.Graph.Inductive.Basic
import Data.Graph.Inductive.Graph
import Data.Graph.Inductive.PatriciaTree
import Data.Graph.Inductive.Query
import Data.List (nub,partition)

import Biobase.Secondary.Diagrams



-- | Given the one to many structures, create the independent graphs, where
-- each graph describes a set of dependent edges in the basepairing.

independentGraphs xs
  | null xs = error "You should give at least one structure. If you did, you have found a bug, yeah!"
  | any ((/=lx).length) xs = error $ "At least one structure with different length was found" ++ show xs
  | otherwise = independentComponents g
  where
    lx = length $ head xs
    g = mkUnionGraph xs

-- | Union of several graphs, created from secondary structure.

mkUnionGraph = unions . map (pairlistToGraph . (mkD1S :: String -> D1Secondary))

-- | Given k graphs, each with nodes [1..n], provide the union of all edges.

unions :: Eq b => [Gr a b] -> Gr a b
unions [] = empty
unions xs@(x:_) = mkGraph (labNodes x) (nub $ concatMap labEdges xs) where

-- | Given a pairlist, generate the secondary structure graph.

pairlistToGraph :: D1Secondary -> Gr () ()
pairlistToGraph d =
  let (len, ps) = fromD1S d
  in  undir $ mkUGraph [0..len -1] ps

-- | Split the graph into (simple, complex) components. The simple components
-- can trivially be filled with any pair. The complex components require an Ear
-- or Woffle decomposition. Simple components are acyclic.

independentComponents :: DynGraph gr => gr () () -> [gr () ()]
independentComponents g = map f $ components g where
  f ns = mkUGraph ns (edges g)

-- | Tests if the given graph is bipartite, which is true if the even/odd BST
-- trees contain no edges

isBipartite :: DynGraph gr => gr a b -> Bool
isBipartite g =  e && o where
  -- True, if there are no edges
  (e,o) = both (null . edges) $ f g
  -- generate the (even,odd) distance graphs
  f g = both (flip delNodes g . map fst) . partition (even . snd) $ level (fst $ head $ labNodes g) g
  both f = second f . first f


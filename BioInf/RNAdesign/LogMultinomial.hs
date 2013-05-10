
module BioInf.RNAdesign.LogMultinomial where

import qualified Data.Vector.Unboxed as VU


logMultinomial :: Int -> VU.Vector Double -> VU.Vector Int -> Double
logMultinomial n' ps xs'
  | VU.length ps /= VU.length xs' = error "logMultinomial: P-vector and count-vector of unequal length"
  | otherwise = n - x + p where
  n = logSum1k n'
  x = VU.sum $ VU.map logSum1k xs'
  p = VU.sum $ VU.zipWith (*) (VU.map log ps) (VU.map fromIntegral xs')
  logSum1k k = VU.sum . VU.map (log . fromIntegral) $ VU.enumFromTo 1 k


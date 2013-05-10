
-- | A simple parser designed to read the optimization string from an argument
-- and together with the input computes the result of the the optimization
-- query. Without 'mkSingleOp' and 'mkMultiOp' this just a trivial parser for
-- simple arithmetic. The addional operations provide access to user-defined
-- functions that can, for example, be used to calculate the energy of a
-- sequence-structure pair. Those functions are not defined here but in the
-- application that uses the parser.

module BioInf.RNAdesign.OptParser
  ( parseOptString
  ) where

import Text.Parsec.Expr
import Text.Parsec hiding ((<|>))
import Text.Parsec.Language
import Text.Parsec.Token
import Control.Applicative
import Text.Parsec.String

import Text.Parsec.Numbers



type SingleOp = (String,Int -> Double)
type MultiOp  = (String,[Double] -> Double)
type GlobalOp = (String,Double)
type PropOp   = (String,[Double] -> Double)
type NumSecStructs = Int

parseOptString :: NumSecStructs -> [SingleOp] -> [MultiOp] -> [GlobalOp] -> [PropOp] -> String -> Double
parseOptString nss sops mops gops props s = case parse expr "" $ prepString s of
  Right res -> res
  Left  err -> error $ show err
  where
    prepString = filter (/= ' ')
    expr :: GenParser Char st Double
    expr
      =   buildExpressionParser optable term
      <?> "expression"
    term :: GenParser Char st Double
    term
      =   between (char '(') (char ')') expr
      <|> parseSingleOp sops
      <|> parseMultiOp nss sops mops
      <|> parseGlobalOp gops
      <|> parsePropOp  props
      <|> parseFloat
      <?> "term"

parsePropOp xs = choice $ map mkPropOp xs

mkPropOp :: PropOp -> GenParser Char st Double
mkPropOp (s,f) = try $ f <$ string s <* string "(" <*> parseFloat `sepBy` string "," <* string ")" where

mkSingleOp :: SingleOp -> GenParser Char st Double
mkSingleOp (s,f) = try $ g <$ string s <* string "(" <*> many1 digit <* string ")" where
  g x = f (read x)

mkMultiOp :: NumSecStructs -> (SingleOp,MultiOp) -> GenParser Char st Double
mkMultiOp nss ((s,sf),(m,mf)) = {- (\xs -> error $ show (xs, map sf xs, mf $ map sf xs)) <$ -} (\xs -> mf $ map sf xs) <$
  string m <* string "(" <* string s <* string "," <*> secs <* string ")" where
    secs  =   try ([1..nss] <$ string "all")
          <|> map read <$> many1 digit `sepBy1` string ","

mkGlobalOp :: GlobalOp -> GenParser Char st Double
mkGlobalOp (s,f) = try $ f <$ string s

parseSingleOp xs = choice $ map mkSingleOp xs

parseMultiOp nss sops mops = choice $ map (try . mkMultiOp nss) [(s,m) | s<-sops, m<-mops]

parseGlobalOp gops = choice $ map (try . mkGlobalOp) gops

optable = [ [prefix "-" negate, prefix "+" id]
          , [binary "^" (**) AssocLeft] --, binary "**" (**) AssocLeft]
          , [binary "*" (*) AssocLeft, binary "/" (/) AssocLeft]
          , [binary "+" (+) AssocLeft, binary "-" (-) AssocLeft]
          ]

pow b e
  | (fromIntegral $ round e) /= e
  = error $ "exponent " ++ show e ++ " needs to be integral, sorry"
  | otherwise
  = b ^ (round e)

prefix name fun       = Prefix (fun <$ string name)
binary name fun assoc = Infix (fun <$ string name) assoc


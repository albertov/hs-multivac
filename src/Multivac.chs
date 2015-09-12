{-# LANGUAGE ForeignFunctionInterface #-}
module Multivac (myMain) where

import Foreign.C.Types (CInt(..))

#include "cbits.h"

{#fun my_main as myMain {} -> `Int' #}

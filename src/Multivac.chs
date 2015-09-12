{-# LANGUAGE ForeignFunctionInterface #-}
module Multivac (
    Mesh
  , newMesh

  , FastMarchSpeedFunc
  , NarrowBandSpeedFunc
  , MaxFSpeedFunc
  , Speed
  , newSpeed

  , simulate
) where

import Control.Applicative ((<$>), (<*>))
import Control.Monad ((>=>))
import Foreign.C.Types
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.ForeignPtr

#include "cbits.h"
{#context prefix = "MV" #}

type FastMarchSpeedFunc  =  Double -> Double -> Double -> IO Double
type NarrowBandSpeedFunc =  Double -> Double -> Double ->
                            Double -> Double -> Double -> IO Double
type MaxFSpeedFunc       =  Double -> Double -> Double ->
                            Double -> Double -> IO Double

{#pointer FastMarchSpeedFunc as CFastMarchSpeedFunc #}
{#pointer NarrowBandSpeedFunc as CNarrowBandSpeedFunc #}
{#pointer MaxFSpeedFunc as CMaxFSpeedFunc #}

foreign import ccall "wrapper"
  wrapFastMarchSpeedFunc :: FastMarchSpeedFunc-> IO CFastMarchSpeedFunc

foreign import ccall "wrapper"
  wrapNarrowBandSpeedFunc :: NarrowBandSpeedFunc-> IO CNarrowBandSpeedFunc

foreign import ccall "wrapper"
  wrapMaxFSpeedFunc :: MaxFSpeedFunc-> IO CMaxFSpeedFunc

{#pointer MeshH as Mesh foreign finalizer DestroyMesh as ^ newtype#}

{#fun NewMesh as ^ {
    `Double'  -- | Xmin
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  } -> `Mesh' #}

{#pointer SpeedH as Speed foreign finalizer DestroySpeed as ^ newtype#}

{#fun NewSpeed as c_NewSpeed {
    `CFastMarchSpeedFunc'
  , `CNarrowBandSpeedFunc'
  , `CMaxFSpeedFunc'
  , `CMaxFSpeedFunc'
  , `Bool'
  , `Bool'
  , `Bool'
  , `Bool'
  } -> `Speed' #}

newSpeed
  :: FastMarchSpeedFunc -> NarrowBandSpeedFunc
  -> MaxFSpeedFunc -> MaxFSpeedFunc -> Bool -> Bool -> Bool -> Bool -> IO Speed
newSpeed fm nb maxF1 maxF2 depPos depTime depNormal depCurv = do
  fm'     <- wrapFastMarchSpeedFunc fm
  nb'     <- wrapNarrowBandSpeedFunc nb
  maxF1'  <- wrapMaxFSpeedFunc maxF1
  maxF2'  <- wrapMaxFSpeedFunc maxF2
  c_NewSpeed fm' nb' maxF1'  maxF2' depPos depTime depNormal depCurv

{#fun Simulate as ^ {`Mesh', `Speed'} -> `()' #}

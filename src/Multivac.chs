{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Multivac (
    Simulation (..)
  , simulate

  , Speed (..)

  , Mesh (..)

  , Front
  , mkFront
  , frontToVector

  , FastMarchSpeedFunc
  , NarrowBandSpeedFunc
  , MaxFSpeedFunc
) where

import Control.Applicative ((<$>), (<*>))
import Control.Monad ((>=>))
import Control.Exception (Exception(toException), SomeException, throw, try)
import Data.Typeable
import Foreign.C.Types
import Foreign.C.String
import Foreign.Marshal.Alloc (alloca)
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import Foreign.ForeignPtr
import System.IO.Unsafe (unsafePerformIO)

import qualified Data.Vector.Storable as St

import Point (Point(..))

data Simulation = Simulation {
    simMesh          :: Mesh
  , simSpeed         :: Speed
  , simNumIterations :: Int
  , simFinalTime     :: Double
}


data Mesh = Mesh {
    meshMinX :: Double
  , meshMinY :: Double
  , meshMaxX :: Double
  , meshMaxY :: Double
  , meshNX   :: Double
  , meshNY   :: Double
  } deriving (Show)

data Speed = Speed {
    speedFastMarch    :: FastMarchSpeedFunc
  , speedNarrowBand   :: NarrowBandSpeedFunc
  , speedMaxF1        :: MaxFSpeedFunc
  , speedMaxF2        :: MaxFSpeedFunc
  , speedDepPosition  :: Bool
  , speedDepTime      :: Bool
  , speedDepNormal    :: Bool
  , speedDepCurvature :: Bool
  }

type ExceptionPtr = StablePtr SomeException

--
-- FastMarchSpeedFunc wrapping
--
type FastMarchSpeedFunc   =  Double -> Double -> Double -> IO Double
type FastMarchSpeedFunc'  =  Double -> Double -> Double -> Ptr CDouble
                          -> IO ExceptionPtr

wrapFastMarchSpeedFunc :: FastMarchSpeedFunc -> IO CFastMarchSpeedFunc
wrapFastMarchSpeedFunc f =
  c_wrapFastMarchSpeedFunc $ \a1 a2 a3 ->
    wrapErrorFunc (f a1 a2 a3) . castPtr

--
-- NarrowBandSpeedFunc wrapping
--
type NarrowBandSpeedFunc  = Double -> Double -> Double -> Double -> Double ->
                            Double -> IO Double
type NarrowBandSpeedFunc' = Double -> Double -> Double -> Double -> Double ->
                            Double -> Ptr CDouble -> IO ExceptionPtr

  
wrapNarrowBandSpeedFunc :: NarrowBandSpeedFunc -> IO CNarrowBandSpeedFunc
wrapNarrowBandSpeedFunc f =
  c_wrapNarrowBandSpeedFunc $ \a1 a2 a3 a4 a5 a6 ->
    wrapErrorFunc (f a1 a2 a3 a4 a5 a6) . castPtr


--
-- MaxSpeedFunc wrapping
--

type MaxFSpeedFunc  =  Double -> Double -> Double -> Double -> Double ->
                       IO Double
type MaxFSpeedFunc' =  Double -> Double -> Double -> Double -> Double ->
                       Ptr CDouble -> IO ExceptionPtr

wrapMaxFSpeedFunc :: MaxFSpeedFunc -> IO CMaxFSpeedFunc
wrapMaxFSpeedFunc f =
  c_wrapMaxFSpeedFunc $ \a1 a2 a3 a4 a5 ->
    wrapErrorFunc (f a1 a2 a3 a4 a5) . castPtr

simulate :: Simulation -> IO ()
simulate Simulation{simMesh=Mesh{..}, simSpeed=Speed{..},..} = do
  speed <- initializeSpeed
  mesh  <- initializeMesh
  c_Simulate mesh speed simNumIterations simFinalTime
  where
    initializeMesh = newMesh meshMinX meshMaxX meshMinY meshMaxY meshNX meshNY
    initializeSpeed = do
      fm    <- wrapFastMarchSpeedFunc speedFastMarch
      nb    <- wrapNarrowBandSpeedFunc speedNarrowBand
      maxF1 <- wrapMaxFSpeedFunc speedMaxF1
      maxF2 <- wrapMaxFSpeedFunc speedMaxF1
      newSpeed fm nb maxF1 maxF2
               speedDepPosition speedDepTime speedDepNormal speedDepCurvature


mkFront :: Orientation -> St.Vector Point -> Front
mkFront orientation v = unsafePerformIO $
  St.unsafeWith v (c_newFront (St.length v) orientation . castPtr)

frontToVector :: Front -> IO (Orientation, St.Vector Point)
frontToVector f = do
  (ptr, size, orientation) <- getFrontPoints f
  fp <- newForeignPtr c_free (castPtr ptr)
  return (orientation, St.unsafeFromForeignPtr0 fp size)



--
-- Error handling
--

data MultivacError = MultivacError String
                   | CallbackError SomeException
  deriving (Show, Typeable)

instance Exception MultivacError

foreign export ccall multivacError :: CString -> IO ExceptionPtr

multivacError :: CString -> IO ExceptionPtr
multivacError = peekCString >=> newStablePtr . toException . MultivacError

wrapErrorFunc
  :: forall a. Storable a
  => IO a -> Ptr a -> IO ExceptionPtr
wrapErrorFunc act ptr = do
    res <- try act
    case res of
      Right v -> poke ptr v >> return nullExceptionPtr
      Left  e -> newStablePtr e
{-# INLINE wrapErrorFunc #-}

nullExceptionPtr :: ExceptionPtr
nullExceptionPtr = castPtrToStablePtr nullPtr

checkHsException :: ExceptionPtr -> IO ()
checkHsException ePtr
  | ePtr == nullExceptionPtr = return ()
  | otherwise                = do exc <- deRefStablePtr ePtr
                                  freeStablePtr ePtr
                                  throw (CallbackError exc)

--
-- Low level bindings
--


#include "cbits.h"
{#context prefix = "MV" #}

{#pointer HsException as ExceptionPtr nocode #}

{#pointer FastMarchSpeedFunc as CFastMarchSpeedFunc #}
{#pointer NarrowBandSpeedFunc as CNarrowBandSpeedFunc #}
{#pointer MaxFSpeedFunc as CMaxFSpeedFunc #}

foreign import ccall "wrapper"
  c_wrapMaxFSpeedFunc :: MaxFSpeedFunc' -> IO CMaxFSpeedFunc

foreign import ccall "wrapper"
  c_wrapFastMarchSpeedFunc :: FastMarchSpeedFunc' -> IO CFastMarchSpeedFunc

foreign import ccall "wrapper"
  c_wrapNarrowBandSpeedFunc :: NarrowBandSpeedFunc' -> IO CNarrowBandSpeedFunc



{#pointer FrontH as Front foreign finalizer DestroyFront as ^ newtype#}

{#enum Orientation {} with prefix = "MVO_" deriving (Eq,Bounded,Show) #}

{#fun NewFront as c_newFront {
    `Int'
  , `Orientation'
  , `Ptr ()'
  } -> `Front' #}

{#fun GetFrontPoints as ^ {
    `Front'
  , alloca- `Int' peekInt*
  , alloca- `Orientation' peekEnum*
  } -> `Ptr ()' #}

peekInt :: Ptr CInt -> IO Int
peekInt = fmap fromIntegral . peek

peekEnum :: Ptr CInt -> IO Orientation
peekEnum = fmap (toEnum . fromIntegral) . peek


{#pointer MeshH foreign finalizer DestroyMesh as ^ newtype#}

{#fun NewMesh as ^ {
    `Double'
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  } -> `MeshH' #}

{#pointer SpeedH foreign finalizer DestroySpeed as ^ newtype#}

{#fun NewSpeed as ^ {
    `CFastMarchSpeedFunc'
  , `CNarrowBandSpeedFunc'
  , `CMaxFSpeedFunc'
  , `CMaxFSpeedFunc'
  , `Bool'
  , `Bool'
  , `Bool'
  , `Bool'
  } -> `SpeedH' #}


{#fun Simulate as c_Simulate {
    `MeshH'
  , `SpeedH'
  , `Int'
  , `Double'
  } -> `()' checkHsException* #}

foreign import ccall "stdlib.h &free"
  c_free :: FunPtr (Ptr a -> IO ())


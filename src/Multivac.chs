{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface #-}
module Multivac (
    Simulation (..)
  , simulate

  , Status (..)

  , Speed (..)

  , Mesh (..)

  , Front
  , mkFront
  , frontToVector
  , throwError

  , FastMarchSpeedFunc
  , NarrowBandSpeedFunc
  , MaxFSpeedFunc
) where

import Control.Applicative ((<$>), (<*>))
import Control.Monad ((>=>), (<=<))
import Control.Exception (Exception(toException), SomeException, throw)
import Data.IORef
import Foreign.C.Types
import Foreign.Marshal.Utils
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

type FastMarchSpeedFunc  =  Double -> Double -> Double -> IO Double
type NarrowBandSpeedFunc =  Double -> Double -> Double ->
                            Double -> Double -> Double -> IO Double
type MaxFSpeedFunc       =  Double -> Double -> Double ->
                            Double -> Double -> IO Double


data Status = Success
            | Error
  deriving (Eq, Show, Bounded, Enum)

simulate :: Simulation -> IO Status
simulate Simulation{simMesh=Mesh{..}, simSpeed=Speed{..},..} = do
  speed <- initializeSpeed
  mesh  <- initializeMesh
  r <- catchAndRethrow (c_Simulate mesh speed simNumIterations simFinalTime)
  return (if r==0 then Success else Error)
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

throwError :: Exception e => e -> IO ()
throwError = (c_throwError . castStablePtrToPtr <=< newStablePtr) . toException

catchAndRethrow :: IO a -> IO a
catchAndRethrow a = do
  result <- newIORef undefined
  ePtr <- c_catchError =<< wrapIOAction (a >>= writeIORef result)
  if ePtr /= nullPtr
    then do
      let sPtr = castPtrToStablePtr ePtr :: StablePtr SomeException
      exc <- deRefStablePtr sPtr
      freeStablePtr sPtr
      throw exc
    else
     readIORef result



--
-- Low level bindings
--


#include "cbits.h"
{#context prefix = "MV" #}

{#pointer FastMarchSpeedFunc as CFastMarchSpeedFunc #}
{#pointer MVAction #}
{#pointer NarrowBandSpeedFunc as CNarrowBandSpeedFunc #}
{#pointer MaxFSpeedFunc as CMaxFSpeedFunc #}

foreign import ccall "wrapper"
  wrapIOAction :: IO () -> IO MVAction

foreign import ccall "wrapper"
  wrapFastMarchSpeedFunc :: FastMarchSpeedFunc-> IO CFastMarchSpeedFunc

foreign import ccall "wrapper"
  wrapNarrowBandSpeedFunc :: NarrowBandSpeedFunc-> IO CNarrowBandSpeedFunc

foreign import ccall "wrapper"
  wrapMaxFSpeedFunc :: MaxFSpeedFunc-> IO CMaxFSpeedFunc

{#fun CatchError as c_catchError {
    `MVAction' }  -> `Ptr ()' #}

{#fun ThrowError as c_throwError {
    `Ptr ()' }  -> `()' #}

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
  } -> `Int' #}

foreign import ccall "stdlib.h &free"
  c_free :: FunPtr (Ptr a -> IO ())


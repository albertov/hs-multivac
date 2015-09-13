{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Multivac (
    Simulation (..)
  , Speed (..)
  , Mesh (..)
  , Front (..)
  , Point (..)
  , Orientation (..)

  , FastMarchSpeedFunc
  , NarrowBandSpeedFunc
  , MaxFSpeedFunc

  , simulate
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
import GHC.Exts (inline)

import qualified Data.Vector.Storable as St
import qualified Data.Vector as V

import Point (Point(..))


data Simulation = Simulation {
    simMesh          :: Mesh
  , simSpeed         :: Speed
  , simNumIterations :: Int
  , simFinalTime     :: Double
  , simInitialFront  :: Front
  , simSavePeriod    :: Double
}

data Front = Front {
    frontCurve       :: St.Vector Point
  , frontOrientation :: Orientation
  } deriving (Show)

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

simulate :: Simulation -> IO (V.Vector Front)
simulate Simulation{..} = do
  fronts <- newFrontArray
  c_Simulate simMesh simSpeed simInitialFront fronts simNumIterations
             simFinalTime simSavePeriod
  toFrontVector fronts
{-# INLINE simulate #-}


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
    wrapErrorFunc (inline f a1 a2 a3) . castPtr
{-# INLINE wrapFastMarchSpeedFunc #-}

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
    wrapErrorFunc (inline f a1 a2 a3 a4 a5 a6) . castPtr
{-# INLINE wrapNarrowBandSpeedFunc #-}


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
    wrapErrorFunc (inline f a1 a2 a3 a4 a5) . castPtr
{-# INLINE wrapMaxFSpeedFunc #-}


toFrontH :: Front -> FrontH
toFrontH f
  | St.null v = error "Empty front"
  | otherwise = unsafePerformIO $ St.unsafeWith v (c_newFront n o . castPtr)
  where Front{frontCurve=v, frontOrientation=o} = f
        n = St.length v

withFront :: Front -> (Ptr FrontH -> IO a) -> IO a
withFront = withFrontH . toFrontH

fromFrontH :: FrontH -> IO Front
fromFrontH f = do
  (ptr, size, orientation) <- getFrontPoints f
  fp <- newForeignPtr c_free (castPtr ptr)
  return $ Front (St.unsafeFromForeignPtr0 fp size) orientation


toFrontVector :: FrontArrayH -> IO (V.Vector Front)
toFrontVector a = do
  n <- getNumFronts a
  front <- c_newEmptyFront
  V.generateM n $ \ix -> do
    copyFrontAt a ix front
    fromFrontH front

toSpeedH :: Speed -> SpeedH
toSpeedH Speed{..} = unsafePerformIO $ do
  fm    <- wrapFastMarchSpeedFunc speedFastMarch
  nb    <- wrapNarrowBandSpeedFunc speedNarrowBand
  maxF1 <- wrapMaxFSpeedFunc speedMaxF1
  maxF2 <- wrapMaxFSpeedFunc speedMaxF1
  newSpeed fm nb maxF1 maxF2
           speedDepPosition speedDepTime speedDepNormal speedDepCurvature

withSpeed :: Speed -> (Ptr SpeedH -> IO a) -> IO a
withSpeed = withSpeedH . toSpeedH

toMeshH :: Mesh -> MeshH
toMeshH Mesh{..} = unsafePerformIO $
  newMesh meshMinX meshMaxX meshMinY meshMaxY meshNX meshNY

withMesh :: Mesh -> (Ptr MeshH -> IO a) -> IO a
withMesh = withMeshH . toMeshH

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


{#pointer FrontH foreign finalizer DestroyFront as ^ newtype#}

{#enum define Orientation { MVO_UNKNOWN       as Unknown
                          , MVO_TRIGONOMETRIC as Trigonometric
                          , MVO_REVERSE       as Reverse
                          } deriving (Eq,Bounded,Show) #}

{#fun unsafe NewFront as c_newFront {
    `Int'
  , `Orientation'
  , `Ptr ()'
  } -> `FrontH' #}

{#fun unsafe NewEmptyFront as c_newEmptyFront  {} -> `FrontH' #}

{#fun unsafe GetFrontPoints as ^ {
    `FrontH'
  , alloca- `Int' peekInt*
  , alloca- `Orientation' peekEnum*
  } -> `Ptr ()' #}

peekInt :: Ptr CInt -> IO Int
peekInt = fmap fromIntegral . peek

peekEnum :: Ptr CInt -> IO Orientation
peekEnum = fmap (toEnum . fromIntegral) . peek


{#pointer FrontArrayH foreign finalizer DestroyFrontArray as ^ newtype#}

{#fun unsafe NewFrontArray as ^ {} -> `FrontArrayH' #}
{#fun unsafe GetNumFronts as ^ {`FrontArrayH'} -> `Int' #}
{#fun unsafe CopyFrontAt as ^ {`FrontArrayH', `Int', `FrontH'} -> `()' #}

{#pointer MeshH foreign finalizer DestroyMesh as ^ newtype#}

{#fun unsafe NewMesh as ^ {
    `Double'
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  , `Double'
  } -> `MeshH' #}

{#pointer SpeedH foreign finalizer DestroySpeed as ^ newtype#}

{#fun unsafe NewSpeed as ^ {
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
    withMesh* `Mesh'
  , withSpeed* `Speed'
  , withFront* `Front'
  , `FrontArrayH'
  , `Int'
  , `Double'
  , `Double'
  } -> `()' checkHsException* #}
{-# INLINE c_Simulate #-}

foreign import ccall "stdlib.h &free"
  c_free :: FunPtr (Ptr a -> IO ())

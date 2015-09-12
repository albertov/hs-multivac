{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface #-}
module Multivac (
    Simulation (..)
  , Status (..)
  , Speed (..)
  , Mesh (..)
  , Front
  , newFront
  , simulate

  , FastMarchSpeedFunc
  , NarrowBandSpeedFunc
  , MaxFSpeedFunc
) where

import Control.Applicative ((<$>), (<*>))
import Control.Monad ((>=>))
import Foreign.C.Types
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.ForeignPtr

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
  r <- c_Simulate mesh speed simNumIterations simFinalTime
  return (if r==0 then Success else Error)
  where
    initializeMesh = newMesh meshMinX meshMinY meshMaxX meshMaxY meshNX meshNY
    initializeSpeed = do
      fm    <- wrapFastMarchSpeedFunc speedFastMarch
      nb    <- wrapNarrowBandSpeedFunc speedNarrowBand
      maxF1 <- wrapMaxFSpeedFunc speedMaxF1
      maxF2 <- wrapMaxFSpeedFunc speedMaxF1
      newSpeed fm nb maxF1 maxF2
               speedDepPosition speedDepTime speedDepNormal speedDepCurvature


newFront :: Orientation -> St.Vector Point -> IO Front
newFront orientation v
  = St.unsafeWith v (c_newFront (St.length v) orientation . castPtr)

--
-- Low level bindings
--


#include "cbits.h"
{#context prefix = "MV" #}

{#pointer FastMarchSpeedFunc as CFastMarchSpeedFunc #}
{#pointer NarrowBandSpeedFunc as CNarrowBandSpeedFunc #}
{#pointer MaxFSpeedFunc as CMaxFSpeedFunc #}

foreign import ccall "wrapper"
  wrapFastMarchSpeedFunc :: FastMarchSpeedFunc-> IO CFastMarchSpeedFunc

foreign import ccall "wrapper"
  wrapNarrowBandSpeedFunc :: NarrowBandSpeedFunc-> IO CNarrowBandSpeedFunc

foreign import ccall "wrapper"
  wrapMaxFSpeedFunc :: MaxFSpeedFunc-> IO CMaxFSpeedFunc

{#pointer FrontH as Front foreign finalizer DestroyFront as ^ newtype#}

{#enum Orientation {} with prefix = "MVO_" deriving (Eq,Bounded,Show) #}

{#fun NewFront as c_newFront {
    `Int'
  , `Orientation'
  , `Ptr ()'
  } -> `Front' #}


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

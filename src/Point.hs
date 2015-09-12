{-# LANGUAGE CPP #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
module Point (Point(..)) where


import Control.Monad (liftM)
import qualified Data.Vector.Generic.Mutable as M
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed.Base as U
import Foreign.Ptr (castPtr)
import Foreign.Storable (Storable(..))

data Point = Point {
    pX :: {-# UNPACK #-} !Double
  , pY :: {-# UNPACK #-} !Double
} deriving (Eq, Show)

instance Storable Point where
  sizeOf _ = 2 * sizeOf (undefined::Double)
  {-# INLINE sizeOf #-}
  alignment _ = alignment (undefined::Double)
  {-# INLINE alignment #-}
  poke ptr (Point x y) = poke ptr' x >> pokeElemOff ptr' 1 y
    where ptr' = castPtr ptr
  {-# INLINE poke #-}
  peek ptr = Point <$> peek ptr' <*> peekElemOff ptr' 1
    where ptr' = castPtr ptr
  {-# INLINE peek #-}

data instance U.Vector    Point =  V_Point !Int (U.Vector    Double)
data instance U.MVector s Point = MV_Point !Int (U.MVector s Double)
instance U.Unbox Point

instance M.MVector U.MVector Point where
  basicLength (MV_Point n _) = n
  basicUnsafeSlice m n (MV_Point _ v) = MV_Point n (M.basicUnsafeSlice (2*m) (2*n) v)
  basicOverlaps (MV_Point _ v) (MV_Point _ u) = M.basicOverlaps v u
  basicUnsafeNew n = liftM (MV_Point n) (M.basicUnsafeNew (2*n))
  basicUnsafeRead (MV_Point _ v) i =
    do let o = 2*i
       x <- M.basicUnsafeRead v o
       y <- M.basicUnsafeRead v (o+1)
       return (Point x y)
  basicUnsafeWrite (MV_Point _ v) i (Point x y) =
    do let o = 2*i
       M.basicUnsafeWrite v o     x
       M.basicUnsafeWrite v (o+1) y
#if MIN_VERSION_vector(0,11,0)
  basicInitialize (MV_Point _ v) = M.basicInitialize v
#endif

instance G.Vector U.Vector Point where
  basicUnsafeFreeze (MV_Point n v) = liftM ( V_Point n) (G.basicUnsafeFreeze v)
  basicUnsafeThaw   ( V_Point n v) = liftM (MV_Point n) (G.basicUnsafeThaw   v)
  basicLength       ( V_Point n _) = n
  basicUnsafeSlice m n (V_Point _ v) = V_Point n (G.basicUnsafeSlice (2*m) (2*n) v)
  basicUnsafeIndexM (V_Point _ v) i =
    do let o = 2*i
       x <- G.basicUnsafeIndexM v o
       y <- G.basicUnsafeIndexM v (o+1)
       return (Point x y)

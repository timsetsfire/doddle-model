package com.picnicml.doddlemodel.robust

package object scale {

  import breeze.numerics._
  import breeze.stats._
  import breeze.stats.distributions.Gaussian

  import com.picnicml.doddlemodel.data.{Features, RealVector, Target}



    object Mad  {

    /** fucntion to calc median absolute deviation
      * median( (X - center(X))/m )
      * http://www.statsmodels.org/dev/_modules/statsmodels/robust/scale.html#mad
      * @param x dense vector of double to be convert to mad
      * @param c normalization constant
      * @param center a funtion used to calculate the c
      */
    def mad(
             x: RealVector,
             c: Double = Gaussian(0,1).inverseCdf(3/4d),
             center: (RealVector => Double) = (x: RealVector) => 0d
    ) = {
      val m = center(x)
      median( abs( x - m) / c)
    }
  }

  case class HuberScale(d: Double = 2.5, tol: Double = 1e-8, maxiter:Double = 30)
}

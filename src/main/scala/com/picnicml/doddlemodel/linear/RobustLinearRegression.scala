package com.picnicml.doddlemodel.linear

import com.picnicml.doddlemodel.data.{Features, RealVector, Target}
import com.picnicml.doddlemodel.robust.norms.{RobustNorm, HuberT}
import breeze.stats.mean

/** An immutable multiple linear regression model with ridge regularization.
  *
  * @param lambda L2 regularization strength, must be positive, 0 means no regularization
  *
  * Examples:
  * val model = LinearRegression()
  * val model = LinearRegression(lambda = 1.5)
  */
@SerialVersionUID(1L)
class RobustLinearRegression private (val rm: RobustNorm = HuberT(), val lambda: Double, protected val w: Option[RealVector])
  extends LinearRegressor[RobustLinearRegression] with Serializable {

  override protected def copy: RobustLinearRegression = new RobustLinearRegression(this.rm, this.lambda, this.w)

  override protected def copy(w: RealVector): RobustLinearRegression = new RobustLinearRegression(this.rm, this.lambda, Some(w))

  @inline override protected def targetVariableAppropriate(y: Target): Boolean = true

  override protected def predict(w: RealVector, x: Features): Target = x * w

  private var yPredCache: Target = _
  private val wSlice: Range.Inclusive = 1 to -1

  override protected[linear] def loss(w: RealVector, x: Features, y: Target): Double = {
    yPredCache = this.predict(w, x)
    val d = rm.rho(y - yPredCache)
    .5 * ((mean(d)) + this.lambda * (w(wSlice).t * w(wSlice)))
  }

  override protected[linear] def lossGrad(w: RealVector, x: Features, y: Target): RealVector = {
    val grad = (rm.psi(y - yPredCache).t * x).t / (-x.rows.toDouble)
    grad(wSlice) += this.lambda * w(wSlice)
    grad
  }
}

object RobustLinearRegression {

  def apply(): RobustLinearRegression = new RobustLinearRegression(lambda = 0, w = None)

  def apply(rm: RobustNorm, lambda: Double): RobustLinearRegression = {
    require(lambda >= 0, "L2 regularization strength must be positive")
    new RobustLinearRegression(rm, lambda, None)
  }
}

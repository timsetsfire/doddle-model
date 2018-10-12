package com.picnicml.doddlemodel.robust

package object norms {

  import breeze.linalg._
  import breeze.numerics._
  import com.picnicml.doddlemodel.data.{Features, RealVector, Target}
  import com.picnicml.doddlemodel.robust.scale._


  //* RobustNorm trait - used for the RlmNet.
  trait RobustNorm {
    def rho(z: RealVector): RealVector
    // the robust ciriterion estiator function
    def psi(z: RealVector): RealVector
    // derivative of rho.  sometimes referred to as the influence function
    def psiDeriv(z: RealVector): RealVector
    // derivative of psi.
    // this can be used to obtain robust covariance matrix
    // for the purpose of enet, will not implement
    def w(z: RealVector): RealVector
    // returns the value of psi(z) / z
  }

  case object L2 extends RobustNorm{
    def rho(z: RealVector) = pow(z, 2) /2d
    def psi(z: RealVector) = z
    def psiDeriv(z: RealVector) = DenseVector.ones[Double](z.length)
    def w(z: RealVector) = DenseVector.ones[Double](z.length)
  }

  case class HuberT(k: Double = 1.345) extends RobustNorm {
    def rho(z: RealVector) = z.map{ elem => if( abs(elem) < k) pow(elem,2) / 2d else k*(abs(elem) - k/2d) }
    def psi(z: RealVector) = z.map{ elem => if( abs(elem) < k) elem else k*signum(elem)}
    def psiDeriv(z: RealVector) = z.map{ elem => if( abs(elem) < k) 1d else 0d}
    def w(z: RealVector) = z.map{ elem => if( abs(elem) < k) 1 else k / abs (elem) }
  }

  case class RamsayE(a: Double = 0.3) extends RobustNorm {
    def rho(z: RealVector): RealVector = {
      val a2inv = pow(a, -2)
      z.map{ elem =>  a2inv * (1 - exp( -a * abs(elem)))*(1d + a * abs(elem)) }
    }
    def psi(z: RealVector): RealVector = z.map{ elem => elem * exp( - a * abs( elem) )}
    def psiDeriv(z: RealVector): RealVector = ???
    def w(z: RealVector): RealVector = exp(abs(z) * (-a))

  }

  case class TukeyBiweight(c: Double  = 4.685) extends RobustNorm {
    def rho(z: RealVector) = z.map{ elem => if(abs(elem) <= c) pow(c,2)/6d * (1d - pow(1d - pow(elem / c,2), 3)) else pow(c,2) / 6d }
    def psi(z: RealVector) = z.map{ elem => if( abs(elem) <= c) elem * pow( 1d - pow(elem/c, 2), 2) else 0d}
    def psiDeriv(z: RealVector) = ???
    def w(z: RealVector) = z.map{ elem => if( abs(elem) <= c) pow( 1d - pow(elem/c, 2), 2) else 0d}
  }

  case class Cauchy(c: Double = 2.3849) extends RobustNorm {
    val c2 = c*c
    def rho(z: RealVector) = z.map{ elem => c2 / 2d * log( 1d + pow(elem/c,2))}
    def psi(z: RealVector) = z.map{ elem => elem / (1d + pow(elem/c,2))}
    def psiDeriv(z: RealVector) = ???
    def w(z: RealVector) = z.map{ elem => 1d / (1d + pow(elem/c,2))}
  }

  case class TrimmedMean(c: Double = 2d) extends RobustNorm {
    def rho(z: RealVector) = z.map{ elem => if( abs(elem) <= c) 1/2d * pow(elem,2) else 0d}
    def psi(z: RealVector) = z.map{ elem => if( abs(elem) <= c) elem else 0d}
    def psiDeriv(z: RealVector) = ???
    def w(z: RealVector) = z.map{ elem => if(abs(elem) <= c) 1 else 0d}
  }

  case object ApproxHuber extends RobustNorm {
    def rho(z: RealVector) = z.map{ elem => 2d* ( sqrt( 1d + pow(elem,2)/2d) - 1d)}
    def psi(z: RealVector) = z.map{ elem => elem / sqrt( 1d + pow(elem,2)/2d) }
    def psiDeriv(z: RealVector) = ???
    def w(z: RealVector) = z.map{ elem => 1d / sqrt( 1d + pow(elem,2)/2d) }
  }

  case object L1 extends RobustNorm {
    def rho(z: RealVector) = abs(z)
    def psi(z: RealVector) = signum(z)
    def psiDeriv(z: RealVector) = ???
    def w(z: RealVector) = z.map{ elem => if(elem == 0d) 1e6 else 1 / abs(elem) }
  }


}

###############################################################################

from pylab import *from numpy import *

###############################################################################

def __unwrap_lms(x, discont=pi, p=10, alpha=0.1, a=[]):
   discont2 = 2*discont
   if a == []:      a = zeros((p))      a[-1] = 1      a = ones((p))/p   for k in xrange(p, x.size):      x_predict = dot(x[k-p:k],a)      e = x[k] - x_predict      e = mod(e+discont, discont2)-discont      x[k] = x_predict + e      a = a + x[k-p:k]*alpha*e   return (x, a)
###############################################################################
def unwrap_lms(x, discont=pi, p=10, alpha=0.1, iter=3):   a = []   start_x = x[0]   for j in xrange(0, 2*iter):      (x, a) = __unwrap_lms(x, discont=discont, p=p, alpha=alpha, a=a)      x = x[::-1].copy()   x = x - x[0] + start_x   return x
###############################################################################
def test_unwrap():   x = matrix(range(0,200)).transpose()/200.0*2*pi   sinus = sin(x)*6   noise = matrix(randn(200,1))*1.1   y = sinus + noise   figure   hold=True;   subplot(5,1,1)   plot(x, sinus, 'k.');   plot(x, y, 'b');   title('Sinusoid in noise');   subplot(5,1,2)   plot(x, mod(sinus,2*pi), 'k.');   plot(x, mod(y,2*pi), 'b');   title('Sinusoid in noise modulo 2pi');   subplot(5,1,3)   plot(x, unwrap(mod(sinus,2*pi),axis=0), 'k.');   plot(x, unwrap(mod(y,2*pi),axis=0), 'b');   title('Unwrapped sinusoid in noise');   subplot(5,1,4)   plot(x, unwrap(mod(sinus,2*pi),axis=0), 'k.');   plot(x, unwrap_lms(y.copy(),p=10,alpha=0,iter=10,axis=0), 'b')   title('Prediction filter unwrapped sinusoid in noise');   subplot(5,1,5)   plot(x, unwrap(mod(sinus,2*pi),axis=0), 'k.');   plot(x, unwrap_lms(y.copy(),p=10,alpha=0.00001,iter=100,axis=0), 'b')   title('Adaptive prediction filter unwrapped sinusoid in noise');if __name__ == "__main__":   test_unwrap()   show()

###############################################################################


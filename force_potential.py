# BERNARD MMAME: STUDENT NUMBER 216075383.
#computational physics project
#finding the force from the potential
import numpy
import random
import matplotlib.pyplot as plt
#
class Potential:
    """ The potential of the sysstem of masses"""
    def __init__(self, m, G, N, soft):
        self.data={}
        self.data['G'] = G
        self.data['N'] = N
        self.data['soft']=soft
        self.data['dt']=dt
        
        self.ax=numpy.zeros(self.data['N'])
        self.ay=numpy.zeros(self.data['N'])
        self.bx=numpy.zeros(self.data['N'])
        self.by=numpy.zeros(self.data['N'])
        self.cx=numpy.zeros(self.data['N'])
        self.cy=numpy.zeros(self.data['N'])
        self.dx=numpy.zeros(self.data['N'])
        self.dy=numpy.zeros(self.data['N'])
        self.lx=[-2,2];self.ly=[0,0]
        self.hx=[0,0];self.hy=[-2,2]
        self.m1=numpy.ones(self.data['N'])*m
        self.m2=self.m1.copy()
        self.m3=self.m1.copy()
        self.m4=self.m1.copy()
        self.M=numpy.sum(self.m1)+numpy.sum(self.m2)+numpy.sum(self.m3)+numpy.sum(self.m4) 
        density = numpy.zeros(self.data['N'])
        #for i in range(len(density)):
         #   density[i]=self.M/16
        
        #self.density = M/16 # 2D which as 4unit dimention 
        for i in range(self.data['N']):
            self.ax[i]=random.uniform(-1.9999,-0.0001)
            self.ay[i]=random.uniform(-1.9999,-0.0001)
            self.bx[i]=random.uniform(0.0001,1.9999)
            self.by[i]=random.uniform(-0.0001,-1.9999)
            self.cx[i]=random.uniform(-0.0001,-1.9999)
            self.cy[i]=random.uniform(0.0001,1.9999)
            self.dx[i]=random.uniform(0.0001,1.9999)
            self.dy[i]=random.uniform(0.0001,1.9999)
            #self.dendity[i] = self.M/16

    
    def Poten_1(self):
        pota=0.0;potb=0.0;potc=0.0;potd=0.0
        for i in range(self.data['N']):
            xa = self.ax[i]-self.ax
            ya = self.ay[i]-self.ay
            ra= xa**2 + ya**2
            soft=self.data['soft']**2
            ra[ra<soft]=soft
            ra=ra+self.data['soft']**2
            Ra = numpy.sqrt(ra)
            pota += self.data['G']*(self.m1/Ra)
            
            xb = 2+self.bx[i]-self.ax
            yb = self.by[i]-self.ay
            rb =xb**2 + yb**2
            soft =self.data['soft']**2
            rb[rb<soft]=soft
            Rb = numpy.sqrt(ra)
            potb  += self.data['G']*(self.m2/Rb)
            
            xc = self.cx[i]-self.ax
            yc =2+self.cy[i]-self.ay
            rc = xc**2 + yc**2
            soft = self.data['soft']**2
            rc[rc<soft]=soft
            Rc = numpy.sqrt(rc)
            potc  += self.data['G']*(self.m3/Rc)

            xd = 2+self.dx[i]-self.ax
            yd = 2+self.dy[i]-self.ay
            rd= xd**2 + yd**2
            soft = self.data['soft']**2
            rd[rd<soft]=soft
            Rd = numpy.sqrt(rd)
            potd  += self.data['G']*(self.m4/Rd)
            return pota+potb+potc+potd
    def Density(self):
        for i in range(self.data['N']):
            density[i] = self.M/16
            #den = density
            return density[i]

    def TotalPoten(self): # is the convolution of dendity and potential 
        d_ft = numpy.fft.fft(sys.Density())
        p_ft = numpy.fft.fft(sys.Poten_1())
        return numpy.real(numpy.fft.ifft(p_ft*d_ft))
    
    def get_force(self): # the gradient of potential
        force = numpy.gradient(sys.TotalPoten())
        return force
        
        
     
        
if __name__=='__main__':
    sys=Potential(1.0,1.0,500,0.01)
    pote=sys.Poten_1()
    Tpote=sys.TotalPoten()
    #print 'the density is: ', den
    print 'the potential of the ith particle is:', pote
    print 'the total potential is ', Tpote
    

plt.plot(sys.ax,sys.ay,'*');plt.plot(sys.bx,sys.by,'*');plt.plot(sys.cx,sys.cy,'*');plt.plot(sys.dx,sys.dy,'*');plt.plot(sys.lx,sys.ly,'r');plt.plot(sys.hx,sys.hy,'r')
plt.show()

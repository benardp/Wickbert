#ifndef KERNEL_H
#define KERNEL_H

#include "Surface/gmTNTConvert.h"
#include "SpaceDef.h"
#include "Surface/Interval.h"
#include "Surface/Box.h"
#include "Surface/KDTree.h"

namespace MLS {

    /**
     * Interface for defining different weight functions for MLS-type surfaces. Kernel 
     * is supposed to be used together with the TGL library, specifically with 
     * ArrayKdTree and ClosePointSet. As such Kernel will store the appropriate
     * ClosePointSet for its own evaluation. For example a finite kernel will store 
     * the kernel size, while an infinite kernel should store the maximal number of
     * neighbors that in practice are used for evaluation
     */
    class Kernel {

    public:
        typedef KDTree<DIMENSION, MLS::Point> KDTree;   ///< Easy to use KDTree defintion

        Kernel() : m_IsFinite(false),m_MaxNumberOfNeighbors(100),m_SqrSize(0),m_Size(0) {}
        virtual ~Kernel() {}

        ///The main interface operation: computing the weight between p1 and p2
        virtual double operator()(const Point &p1, const Point &p2) = 0;

        ///For interval analysis we may want to compute the ranges of weights between a point
        ///and a box.
        virtual Interval<double> operator() (const Box<double> &b, Point &p) = 0;

        ///First order partial differential equation of the weight over the box b to point
        ///p.  Differentiating with respect to arg.
        virtual Interval<double> pde1Box(const Box<double> &b, Point &p, int arg) = 0;

        bool isFinite() {return m_IsFinite;}
		void setIsFinite(bool b) {m_IsFinite = b;}
        int getMaxNeighbors() {return m_MaxNumberOfNeighbors;}
        double getKernelSizeSquared() {return m_SqrSize;}
        double getMahFactor() {return m_MahFactor;}
        void setMaxNeighbors(int neigh) {m_MaxNumberOfNeighbors = neigh;}
        void setKernelSize(double size) {m_SqrSize = size;}
        void setMahFactor(double f) {m_MahFactor = f;}
    
        /**
         * Gets the neighbors of this point as the Kernel sees fit.
         */
        virtual int getNeighbors(const Point &query, KDTree *kdTree,
            std::vector<DistancedIndex> &neighbors) {

            if (this->isFinite())
                return kdTree->nearestInRadiusQuery(query, m_Size, m_MaxNumberOfNeighbors, neighbors);
            else
                return kdTree->nearestQuery(query, m_MaxNumberOfNeighbors, neighbors);
        }

        /**
         * Gets the neighbors of this point as the Kernel sees fit.
         */

        virtual int getNeighbors(const Box<double> &query, KDTree *kdTree,
            std::vector<DistancedIndex> &neighbors) {

            //Find an appropriate radius to query.  This is the interval diagonal/2 + kernel finite radius.
            if (this->isFinite())
                return kdTree->nearestInRadiusQuery(convert(query.center()), 
                    query.diagonal()/2. + m_Size, m_MaxNumberOfNeighbors, neighbors);
            else
                return kdTree->nearestQuery(convert(query.center()), m_MaxNumberOfNeighbors, neighbors);
        }

    protected:
        ///Constructors only to be used by derived classes
        Kernel(int neighSize) : m_IsFinite(false),m_MaxNumberOfNeighbors(neighSize),m_SqrSize(0),m_Size(0) {}
        Kernel(int neighSize, double size) : m_IsFinite(true),m_MaxNumberOfNeighbors(neighSize),m_SqrSize(size) {
            m_Size = sqrt(m_SqrSize);
        }
        
        /**
         * Infinite Kernel, size here is used to indicate how far from 0 the kernel is non-zero
         * (sqrSize is this squared). HACK ALERT.  0 is the sentinel value for KDTree to give you
         * no limit to the number of entries.
         */
        Kernel(double sqrSize) : m_IsFinite(true),m_MaxNumberOfNeighbors(0),m_SqrSize(sqrSize) {}

    protected:

        bool m_IsFinite;            ///< is the kernel finite or infinite
        int m_MaxNumberOfNeighbors; ///< the maximal number of neighbors used for evaluation
        double m_SqrSize;           ///< the squared size of a finite kernel
        double m_Size;              ///< Sqrt of m_SqrSize
        double m_MahFactor;
    };

    class GaussKernel : public Kernel {

    public:
        typedef Kernel BaseKernel;
        
        GaussKernel(int neighSize = 15, double var = 0.25) : BaseKernel(neighSize) {m_SqrVariance = var;}
        GaussKernel(int neighSize, double radius2, double var) : BaseKernel(neighSize, radius2) {m_SqrVariance = var;}
        ~GaussKernel() {}

        ///The main interface operation: computing the weight between p1 and p2
        virtual double operator() (const Point &p1, const Point &p2) {
            
            double dist2 = distanceSquared(p1,p2);
            return (isFinite() && (dist2 > m_SqrSize)) ? 0 : exp((-dist2) / m_SqrVariance);
        }

       ///For interval analysis we may want to compute the ranges of weights between a point
       ///and a box.
		
       virtual Interval<double> operator() (const Box<double> &b, Point &p) {

            Interval<double> top = (b[0] - p[0]).squared() + (b[1] - p[1]).squared() + (b[2] - p[2]).squared();
            return (top/(-m_SqrVariance)).exp();
        }
		
        ///First order partial differential equation of the weight over the box b to point
        ///p.  Differentiating with respect to arg.
		
       virtual Interval<double> pde1Box(const Box<double> &b, Point &p, int arg) {

            Interval<double> top = (b[0] - p[0]).squared() + (b[1] - p[1]).squared() + (b[2] - p[2]).squared();
            top /= -m_SqrVariance;
            return (-2/m_SqrVariance)*(b[arg] - p[arg]) * top.exp();
       }
	   

		double getSqrVariance() {return m_SqrVariance;}
        void setSqrVariance(double s) {m_SqrVariance = s;}

    protected:

        double m_SqrVariance;
    };

    /**
     * Kernel with a Gaussian center and linear tails to 0.
     */
    class GaussianClip : public Kernel {

    public:
        typedef Kernel BaseKernel;
        
        GaussianClip(double radius2 = 0.25, double var2 = 0.25, double linearCO = 0.1) : BaseKernel(radius2) {
            m_SqrVariance = var2;
            m_Size = sqrt(radius2);

            m_LinearCutOff = m_Size * linearCO;
            m_LinearSlope = -1 / (m_Size - m_LinearCutOff);
            m_LinearOffset = -(m_LinearSlope*m_Size);
        }
        ~GaussianClip() {}

        ///The main interface operation: computing the weight between p1 and p2
        virtual double operator() (const Point &p1, const Point &p2) {
            
            double dist = distance(p1, p2);
            assert(dist <= m_Size);

            if (dist < m_LinearCutOff) {
                return exp(-(dist*dist)) / m_SqrVariance;
            }
            else {
                return (exp(-(dist*dist)) / m_SqrVariance) * (dist*m_LinearSlope + m_LinearOffset);
            }
        }

       ///For interval analysis we may want to compute the ranges of weights between a point
       ///and a box.
		
       virtual Interval<double> operator() (const Box<double> &b, Point &p) {

            Interval<double> top = (b[0] - p[0]).squared() + (b[1] - p[1]).squared() + (b[2] - p[2]).squared();
            Interval<double> retVal = (top/(-m_SqrVariance)).exp();

            //Now we need to attenuate the low end of the interval even further if it is
            //in the tail of the gaussian that gets clipped
            double bigD = sqrt(top.high());
            double smallD = sqrt(top.low());
            double high = (smallD > m_LinearCutOff) ? 
                MAX(retVal.high() * (smallD * m_LinearSlope + m_LinearOffset), 0.) : retVal.high();
            double low = (bigD > m_LinearCutOff) ? 
                MAX(retVal.low() * (bigD * m_LinearSlope + m_LinearOffset), 0.) : retVal.low();
            retVal.setInterval(low, high);

            return retVal;
       }
		
        ///First order partial differential equation of the weight over the box b to point
        ///p.  Differentiating with respect to arg.
		
        virtual Interval<double> pde1Box(const Box<double> &b, Point &p, int arg) {
           
            Interval<double> dist2 = (b[0] - p[0]).squared() + (b[1] - p[1]).squared() + (b[2] - p[2]).squared();
            Interval<double> top = dist2 / -m_SqrVariance;
            top /= -m_SqrVariance;
            Interval<double> argInt = b[arg] - p[arg];
            Interval<double> retVal =(-2/m_SqrVariance)*argInt * top.exp();

            //We are only worried about cases where the distance to the point makes us use our clipping
            //piece.  In this case the derivative might have greater magnitude than the above calculation
            //indicates.  The magnitude will never be too small so long as argInt contains 0, because
            //then we're getting the derivative of the steepest part of the Gaussian.  All we handle
            //below are cases where we are far from the point and we don't necessarily contain the steepest
            //part of the gaussian.
            if ((dist2.high() > (m_LinearCutOff * m_LinearCutOff)) && (!argInt.contains(0))) {
             
                if(argInt.isPositive()) {

                    double newLow = \
						(m_LinearSlope*exp((m_LinearCutOff*m_LinearCutOff) / -m_SqrVariance)) \
						+ \
						(exp((m_LinearCutOff*m_LinearCutOff) / -m_SqrVariance)*(-2/m_SqrVariance)*argInt.high()) \
						* \
                            (m_LinearCutOff*m_LinearSlope + m_LinearOffset);
                    //only set the lower part of the derivative interval, as it will have greater magnitude
                    retVal.setInterval(newLow, retVal.high());
                }
                else {

                    double newHigh = \
                        (m_LinearSlope*exp((m_LinearCutOff*m_LinearCutOff) / -m_SqrVariance)) \
                        + \
						(exp((m_LinearCutOff*m_LinearCutOff) / -m_SqrVariance)*(-2/m_SqrVariance)*argInt.low()) \
						* \
						(m_LinearCutOff*m_LinearSlope + m_LinearOffset);
                    //only set the higher part ofthe derivative interval, as it will have greater magnitude
                    retVal.setInterval(retVal.low(), newHigh);
                }
            }

            return retVal;
        }
		
        /**
         * Gets the neighbors of this point as the Kernel sees fit.
         */
        virtual int getNeighbors(const Point &query, KDTree *kdTree,
            std::vector<DistancedIndex> &neighbors) {

            return kdTree->radiusQuery(query, m_Size, neighbors);
        }

        /**
         * Gets the neighbors of this point as the Kernel sees fit.
         */
		
        virtual int getNeighbors(const Box<double> &query, KDTree *kdTree,
            std::vector<DistancedIndex> &neighbors) {

            //Find an appropriate radius to query.  This is the interval diagonal/2 + kernel finite radius.
            return kdTree->radiusQuery(convert(query.center()),
                m_Size + query.diagonal()/2, neighbors);
        }
	
    public:
		double getSqrVariance() {return m_SqrVariance;}
        void setSqrVariance(double s) {m_SqrVariance = s;}
		double getLinearCutOff() {return m_LinearCutOff;}
        void setLinearCutOff(double s) {m_LinearCutOff = s;}

    private:
		double m_SqrVariance;   ///< Square of the variance in the Gaussian.
        double m_LinearCutOff;  ///< When the linear clipper kicks in (asbolute distance from 0)
        double m_LinearOffset;  ///< y = mx + B
        double m_LinearSlope;   ///< y = Mx + b;
    };
};

#endif // KERNEL_H

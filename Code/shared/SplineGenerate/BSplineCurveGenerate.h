// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Mathematics/BasisFunction.h>
#include <Mathematics/BandedMatrix.h>
#include <iostream>
using namespace std;

// The algorithm implemented here is based on the document
// https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf

namespace gte
{
    template <typename Real>
    class BSplineCurveGenerate
    {
    public:
        // Construction.  The preconditions for calling the constructor are
        //   1 <= degree && degree < numControls <= numSamples
        // The samples points are contiguous blocks of 'dimension' real values
        // stored in sampleData.
        BSplineCurveGenerate(int dimension, int degree, vector<Real> ControlData, int numControls)
            :
            mDimension(dimension),
            mDegree(degree),
            mControlData(ControlData),
            mNumControls(numControls)
        {
            LogAssert(dimension >= 1, "Invalid dimension.");
            LogAssert(1 <= degree && degree < numControls, "Invalid degree.");
           
            BasisFunctionInput<Real> input;
            input.numControls = numControls;
            input.degree = degree;
            input.uniform = true;
            input.periodic = false;
            input.numUniqueKnots = numControls - degree + 1;
            input.uniqueKnots.resize(input.numUniqueKnots);
            input.uniqueKnots[0].t = (Real)0;
            input.uniqueKnots[0].multiplicity = degree + 1;
            int last = input.numUniqueKnots - 1;
            Real factor = ((Real)1) / (Real)last;
            for (int i = 1; i < last; ++i)
            {
                input.uniqueKnots[i].t = factor * (Real)i;
                input.uniqueKnots[i].multiplicity = 1;
            }
            input.uniqueKnots[last].t = (Real)1;
            input.uniqueKnots[last].multiplicity = degree + 1;
            mBasis.Create(input);
        }

        // Access to input sample information.
        inline int GetDimension() const
        {
            return mDimension;
        }

        inline int GetNumSamples() const
        {
            return mNumSamples;
        }

        inline Real const* GetSampleData() const
        {
            return mSampleData;
        }

        // Access to output control point and curve information.
        inline int GetDegree() const
        {
            return mDegree;
        }

        inline int GetNumControls() const
        {
            return mNumControls;
        }

        inline Real const* GetControlData() const
        {
            return &mControlData[0];
        }

        inline BasisFunction<Real> const& GetBasis() const
        {
            return mBasis;
        }

        // Evaluation of the B-spline curve.  It is defined for 0 <= t <= 1.
        // If a t-value is outside [0,1], an open spline clamps it to [0,1].
        // The caller must ensure that position[] has at least 'dimension'
        // elements.
        void Evaluate(Real t, unsigned int order, Real* value) const
        {
            int imin, imax;
            mBasis.Evaluate(t, order, imin, imax);

            Real const* source = &mControlData[mDimension * imin];
            Real basisValue = mBasis.GetValue(order, imin);
            for (int j = 0; j < mDimension; ++j)
            {
                value[j] = basisValue * (*source++);
            }

            for (int i = imin + 1; i <= imax; ++i)
            {
                basisValue = mBasis.GetValue(order, i);
                for (int j = 0; j < mDimension; ++j)
                {
                    value[j] += basisValue * (*source++);
                }
            }
        }

        void GetPosition(Real t, Real* position) const
        {
            Evaluate(t, 0, position);
        }
/*
        void SplineGenerate(Real t, vector<Real> ControlData, Real* position) const
        {
            int imin, imax;
            mBasis.Evaluate(t, 0, imin, imax);

            Real const* source = &ControlData[mDimension * imin];
            Real basisValue = mBasis.GetValue(0, imin);
            for (int j = 0; j < mDimension; ++j)
            {
                position[j] = basisValue * (*source++);
            }

            for (int i = imin + 1; i <= imax; ++i)
            {
                basisValue = mBasis.GetValue(0, i);
                for (int j = 0; j < mDimension; ++j)
                {
                    position[j] += basisValue * (*source++);
                }
            }
        }
        */
    private:
        // Input sample information.
        int mDimension;
        int mNumSamples;
        Real const* mSampleData;

        // The fitted B-spline curve, open and with uniform knots.
        int mDegree;
        std::vector<Real> mControlData;
        int mNumControls;
        BasisFunction<Real> mBasis;
    };
}

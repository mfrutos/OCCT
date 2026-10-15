// Created on: 1995-07-18
// Created by: Modelistation
// Copyright (c) 1995-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#ifndef _Extrema_GGenExtCC_HeaderFile
#define _Extrema_GGenExtCC_HeaderFile

#include <algorithm>

#include <Extrema_GlobOptFuncCC.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <math_GlobOptMin.hxx>
#include <math_Vector.hxx>
#include <NCollection_CellFilter.hxx>
#include <NCollection_DynamicArray.hxx>
#include <Precision.hxx>
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_NullObject.hxx>
#include <Standard_OutOfRange.hxx>
#include <StdFail_NotDone.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <Standard_Integer.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Sequence.hxx>
#include <ElCLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <gp_Circ.hxx>
#include <gp_Circ2d.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <gp_Vec2d.hxx>
#include <gp_XYZ.hxx>

#include <cmath>
#include <functional>
#include <queue>
#include <vector>

#ifndef M_SQRT2
  #define M_SQRT2 1.41421356237309504880168872420969808
#endif

//! Template class for computing extremal distances between two curves.
//! The function F(u,v)=distance(C1(u),C2(v)) has an extremum when gradient(f)=0.
//! The algorithm uses Evtushenko's global optimization solver.
//!
//! @tparam TheCurve1 Type of the first curve (e.g., Adaptor3d_Curve)
//! @tparam TheCurveTool1 Tool class for the first curve
//! @tparam TheCurve2 Type of the second curve
//! @tparam TheCurveTool2 Tool class for the second curve
//! @tparam ThePOnC Point on curve type (e.g., Extrema_POnCurv)
//! @tparam ThePoint Point type (e.g., gp_Pnt)
//! @tparam TheExtPC Point-to-curve extremum class (e.g., Extrema_ExtPC)
template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
class Extrema_GGenExtCC
{
public:
  DEFINE_STANDARD_ALLOC

  //! Default constructor.
  Extrema_GGenExtCC();

  //! Constructor with two curves.
  Extrema_GGenExtCC(const TheCurve1& theC1, const TheCurve2& theC2);

  //! Constructor with two curves and parameter bounds.
  Extrema_GGenExtCC(const TheCurve1& theC1,
                    const TheCurve2& theC2,
                    const double     theUinf,
                    const double     theUsup,
                    const double     theVinf,
                    const double     theVsup);

  //! Sets parameters for computation.
  void SetParams(const TheCurve1& theC1,
                 const TheCurve2& theC2,
                 const double     theUinf,
                 const double     theUsup,
                 const double     theVinf,
                 const double     theVsup);

  //! Sets the tolerance.
  void SetTolerance(const double theTol);

  //! Set flag for single extrema computation.
  void SetSingleSolutionFlag(const bool theFlag);

  //! Get flag for single extrema computation.
  bool GetSingleSolutionFlag() const;

  //! Performs calculations.
  void Perform();

  //! Returns True if the distances are found.
  bool IsDone() const;

  //! Returns state of myParallel flag.
  bool IsParallel() const;

  //! Returns the number of extremum distances.
  int NbExt() const;

  //! Returns the value of the Nth square extremum distance.
  double SquareDistance(const int theN = 1) const;

  //! Returns the points of the Nth extremum distance.
  void Points(const int theN, ThePOnC& theP1, ThePOnC& theP2) const;

private:
  bool                         myIsFindSingleSolution;
  bool                         myParallel;
  double                       myCurveMinTol;
  math_Vector                  myLowBorder;
  math_Vector                  myUppBorder;
  NCollection_Sequence<double> myPoints1;
  NCollection_Sequence<double> myPoints2;
  void*                        myC[2];
  bool                         myDone;
};

//=================================================================================================
// Static helper functions
//=================================================================================================

namespace
{
// Comparator, used in std::sort.
inline bool Extrema_GGenExtCC_comp(const gp_XY& theA, const gp_XY& theB)
{
  if (theA.X() < theB.X())
  {
    return true;
  }
  else
  {
    if (theA.X() == theB.X())
    {
      if (theA.Y() < theB.Y())
        return true;
    }
  }
  return false;
}

inline void Extrema_GGenExtCC_ChangeIntervals(occ::handle<NCollection_HArray1<double>>& theInts,
                                              const int                                 theNbInts)
{
  int                                      aNbInts = theInts->Length() - 1;
  int                                      aNbAdd  = theNbInts - aNbInts;
  occ::handle<NCollection_HArray1<double>> aNewInts =
    new NCollection_HArray1<double>(1, theNbInts + 1);
  int aNbLast = theInts->Length();
  int i;
  if (aNbInts == 1)
  {
    aNewInts->SetValue(1, theInts->First());
    aNewInts->SetValue(theNbInts + 1, theInts->Last());
    double dt = (theInts->Last() - theInts->First()) / theNbInts;
    double t  = theInts->First() + dt;
    for (i = 2; i <= theNbInts; ++i, t += dt)
    {
      aNewInts->SetValue(i, t);
    }
    theInts = aNewInts;
    return;
  }
  for (i = 1; i <= aNbLast; ++i)
  {
    aNewInts->SetValue(i, theInts->Value(i));
  }
  while (aNbAdd > 0)
  {
    double anLIntMax = -1.;
    int    aMaxInd   = -1;
    for (i = 1; i < aNbLast; ++i)
    {
      double anL = aNewInts->Value(i + 1) - aNewInts->Value(i);
      if (anL > anLIntMax)
      {
        anLIntMax = anL;
        aMaxInd   = i;
      }
    }

    double t = (aNewInts->Value(aMaxInd + 1) + aNewInts->Value(aMaxInd)) / 2.;
    for (i = aNbLast; i > aMaxInd; --i)
    {
      aNewInts->SetValue(i + 1, aNewInts->Value(i));
    }
    aNbLast++;
    aNbAdd--;
    aNewInts->SetValue(aMaxInd + 1, t);
  }
  theInts = aNewInts;
}

class Extrema_GGenExtCC_PointsInspector
{
public:
  static constexpr int Dimension = 2;

  typedef gp_XY Point;
  typedef gp_XY Target;

  static double Coord(int i, const Point& thePnt) { return thePnt.Coord(i + 1); }

  static Point Shift(const Point& thePnt, double theTol)
  {
    return Point(thePnt.X() + theTol, thePnt.Y() + theTol);
  }

  Extrema_GGenExtCC_PointsInspector(const double theTol)
      : myTol(theTol * theTol),
        myIsFind(false)
  {
  }

  void ClearFind() { myIsFind = false; }

  bool isFind() const { return myIsFind; }

  void SetCurrent(const gp_XY& theCurPnt) { myCurrent = theCurPnt; }

  NCollection_CellFilter_Action Inspect(const Target& theObject)
  {
    gp_XY        aPt     = myCurrent.Subtracted(theObject);
    const double aSQDist = aPt.SquareModulus();
    if (aSQDist < myTol)
    {
      myIsFind = true;
    }
    return CellFilter_Keep;
  }

private:
  double myTol;
  gp_XY  myCurrent;
  bool   myIsFind;
};

template <typename TheCurve, typename TheExtPCType, typename ThePointType>
double Extrema_GGenExtCC_ProjPOnC(const ThePointType& theP, TheExtPCType& theProjTool)
{
  double aDist = ::RealLast();
  theProjTool.Perform(theP);
  if (theProjTool.IsDone() && theProjTool.NbExt())
  {
    for (int i = 1; i <= theProjTool.NbExt(); ++i)
    {
      double aD = theProjTool.SquareDistance(i);
      if (aD < aDist)
        aDist = aD;
    }
  }
  return aDist;
}

// A branch and bound that certifies the least distance between lines, circles,
// Bezier and B-spline curves from exact bounds on |C'| and |C''|.
struct Extrema_CertCC_Piece
{
  double Lo, Hi, L1, L2;
};

struct Extrema_CertCC_Homog
{
  gp_XYZ P;
  double W;
};

inline gp_XYZ Extrema_CertCC_XYZ(const gp_Pnt& theP)
{
  return theP.XYZ();
}

inline gp_XYZ Extrema_CertCC_XYZ(const gp_Pnt2d& theP)
{
  return gp_XYZ(theP.X(), theP.Y(), 0.0);
}

template <typename ThePoint>
struct Extrema_CertCC_VecOf;

template <>
struct Extrema_CertCC_VecOf<gp_Pnt>
{
  using Type = gp_Vec;
};

template <>
struct Extrema_CertCC_VecOf<gp_Pnt2d>
{
  using Type = gp_Vec2d;
};

// One derivative of the homogeneous control points theC of one span, of
// degree theP, whose first control point has index theI0 in the flat knots theT.
inline std::vector<Extrema_CertCC_Homog> Extrema_CertCC_Derive(
  const std::vector<Extrema_CertCC_Homog>& theC,
  const int                                theP,
  const int                                theI0,
  const std::vector<double>&               theT,
  const int                                theLevel)
{
  std::vector<Extrema_CertCC_Homog> aQ;
  for (size_t j = 0; j + 1 < theC.size(); ++j)
  {
    const int    i    = theI0 + static_cast<int>(j);
    const double aDen = theT[i + theP + theLevel] - theT[i + theLevel];
    if (aDen > 0.0)
    {
      const double aS = theP / aDen;
      aQ.push_back({(theC[j + 1].P - theC[j].P) * aS, (theC[j + 1].W - theC[j].W) * aS});
    }
    else
    {
      aQ.push_back({gp_XYZ(0.0, 0.0, 0.0), 0.0});
    }
  }
  return aQ;
}

// Cuts a non-periodic spline (poles, weights, flat knots) into its spans over
// [theLo, theHi], repeating the spans every thePeriod when it is positive.
// For a rational span, with H = w (P - O): C' = (H' - (C - O) w') / w and
// C'' = (H'' - 2 C' w' - (C - O) w'') / w, and |C - O| <= max |P_i - O|.
inline bool Extrema_CertCC_SplinePieces(const std::vector<gp_XYZ>&         thePoles,
                                        const std::vector<double>&         theWeights,
                                        const std::vector<double>&         theT,
                                        const int                          theDeg,
                                        const double                       thePeriod,
                                        const double                       theLo,
                                        const double                       theHi,
                                        std::vector<Extrema_CertCC_Piece>& thePieces)
{
  const int p = theDeg;
  const int n = static_cast<int>(thePoles.size());
  if (p < 1 || static_cast<int>(theT.size()) != n + p + 1)
  {
    return false;
  }
  int aLowShift = 0, aHighShift = 0;
  if (thePeriod > 0.0)
  {
    aLowShift  = static_cast<int>(std::floor((theLo - theT[p]) / thePeriod)) - 1;
    aHighShift = static_cast<int>(std::ceil((theHi - theT[p]) / thePeriod)) + 1;
  }
  for (int k = p; k < n; ++k)
  {
    if (!(theT[k + 1] > theT[k]))
    {
      continue;
    }
    const gp_XYZ&                     anO = thePoles[k - p];
    std::vector<Extrema_CertCC_Homog> aH;
    double                            aWMin = RealLast(), aR = 0.0;
    for (int i = k - p; i <= k; ++i)
    {
      const double aW = theWeights[i];
      aH.push_back({(thePoles[i] - anO) * aW, aW});
      aWMin = std::min(aWMin, aW);
      aR    = std::max(aR, (thePoles[i] - anO).Modulus());
    }
    if (!(aWMin > 0.0))
    {
      return false;
    }
    const std::vector<Extrema_CertCC_Homog> aH1 = Extrema_CertCC_Derive(aH, p, k - p, theT, 1);
    const std::vector<Extrema_CertCC_Homog> aH2 =
      p >= 2 ? Extrema_CertCC_Derive(aH1, p - 1, k - p, theT, 2)
             : std::vector<Extrema_CertCC_Homog>();
    double anA1 = 0.0, aW1 = 0.0, anA2 = 0.0, aW2 = 0.0;
    for (const Extrema_CertCC_Homog& aQ : aH1)
    {
      anA1 = std::max(anA1, aQ.P.Modulus());
      aW1  = std::max(aW1, std::abs(aQ.W));
    }
    for (const Extrema_CertCC_Homog& aQ : aH2)
    {
      anA2 = std::max(anA2, aQ.P.Modulus());
      aW2  = std::max(aW2, std::abs(aQ.W));
    }
    const double aL1 = (anA1 + aR * aW1) / aWMin;
    const double aL2 = (anA2 + 2.0 * aL1 * aW1 + aR * aW2) / aWMin;
    for (int aShift = aLowShift; aShift <= aHighShift; ++aShift)
    {
      const double aS0 = std::max(theT[k] + aShift * thePeriod, theLo);
      const double aS1 = std::min(theT[k + 1] + aShift * thePeriod, theHi);
      if (aS1 > aS0)
      {
        thePieces.push_back({aS0, aS1, aL1, aL2});
      }
    }
  }
  std::sort(thePieces.begin(),
            thePieces.end(),
            [](const Extrema_CertCC_Piece& theA, const Extrema_CertCC_Piece& theB) {
              return theA.Lo < theB.Lo;
            });
  // A span shifted by whole periods ends within rounding of the next one.
  for (size_t i = 1; i < thePieces.size(); ++i)
  {
    if (std::abs(thePieces[i].Lo - thePieces[i - 1].Hi)
        <= 1.0e-12 * (1.0 + std::abs(thePieces[i].Lo)))
    {
      thePieces[i].Lo = thePieces[i - 1].Hi;
    }
  }
  // Every cell is cut at piece ends, so a gap would leave part of the range unbounded.
  bool isCovered =
    !thePieces.empty() && thePieces.front().Lo == theLo && thePieces.back().Hi == theHi;
  for (size_t i = 1; isCovered && i < thePieces.size(); ++i)
  {
    isCovered = thePieces[i].Lo == thePieces[i - 1].Hi;
  }
  return isCovered;
}

// The bounds are taken from the geometry the adaptor reports. An adaptor whose
// BSpline() is an approximation (a curve on a surface) does not evaluate to it.
template <typename TheCurve, typename TheCurveTool, typename TheGeom>
bool Extrema_CertCC_SameCurve(const TheCurve& theC,
                              const TheGeom&  theGeom,
                              const double    theLo,
                              const double    theHi)
{
  for (const double aT : {theLo, 0.5 * (theLo + theHi), theHi})
  {
    const gp_XYZ aP = Extrema_CertCC_XYZ(theGeom(aT));
    const gp_XYZ aQ = Extrema_CertCC_XYZ(TheCurveTool::Value(theC, aT));
    if ((aP - aQ).Modulus() > 1.0e-12 * (1.0 + aP.Modulus()))
    {
      return false;
    }
  }
  return true;
}

template <typename TheCurve, typename TheCurveTool>
bool Extrema_CertCC_Pieces(const TheCurve&                    theC,
                           const double                       theLo,
                           const double                       theHi,
                           std::vector<Extrema_CertCC_Piece>& thePieces)
{
  if (!(theHi > theLo) || Precision::IsInfinite(theLo) || Precision::IsInfinite(theHi))
  {
    return false;
  }
  switch (TheCurveTool::GetType(theC))
  {
    case GeomAbs_Line: {
      const auto aLin = TheCurveTool::Line(theC);
      if (!Extrema_CertCC_SameCurve<TheCurve, TheCurveTool>(
            theC,
            [&](const double theT) { return ElCLib::Value(theT, aLin); },
            theLo,
            theHi))
      {
        return false;
      }
      thePieces.push_back({theLo, theHi, 1.0, 0.0});
      return true;
    }
    case GeomAbs_Circle: {
      const auto aCirc = TheCurveTool::Circle(theC);
      if (!Extrema_CertCC_SameCurve<TheCurve, TheCurveTool>(
            theC,
            [&](const double theT) { return ElCLib::Value(theT, aCirc); },
            theLo,
            theHi))
      {
        return false;
      }
      const double aR  = aCirc.Radius();
      const int    aNb = std::max(1, static_cast<int>(std::ceil((theHi - theLo) / (M_PI / 4.0))));
      for (int i = 0; i < aNb; ++i)
      {
        thePieces.push_back({theLo + (theHi - theLo) * i / aNb,
                             i + 1 == aNb ? theHi : theLo + (theHi - theLo) * (i + 1) / aNb,
                             aR,
                             aR});
      }
      return true;
    }
    case GeomAbs_BezierCurve: {
      const auto aBez = TheCurveTool::Bezier(theC);
      if (aBez.IsNull()
          || !Extrema_CertCC_SameCurve<TheCurve, TheCurveTool>(
            theC,
            [&](const double theT) { return aBez->Value(theT); },
            theLo,
            theHi))
      {
        return false;
      }
      const int           aDeg = aBez->Degree();
      std::vector<gp_XYZ> aPoles;
      std::vector<double> aWeights, aT;
      for (int i = 1; i <= aBez->NbPoles(); ++i)
      {
        aPoles.push_back(Extrema_CertCC_XYZ(aBez->Pole(i)));
        aWeights.push_back(aBez->Weight(i));
      }
      aT.assign(aDeg + 1, 0.0);
      aT.insert(aT.end(), aDeg + 1, 1.0);
      return Extrema_CertCC_SplinePieces(aPoles, aWeights, aT, aDeg, 0.0, theLo, theHi, thePieces);
    }
    case GeomAbs_BSplineCurve: {
      const auto aBSpl = TheCurveTool::BSpline(theC);
      if (aBSpl.IsNull()
          || !Extrema_CertCC_SameCurve<TheCurve, TheCurveTool>(
            theC,
            [&](const double theT) { return aBSpl->Value(theT); },
            theLo,
            theHi))
      {
        return false;
      }
      auto   aFlat   = aBSpl;
      double aPeriod = 0.0;
      if (aBSpl->IsPeriodic())
      {
        aPeriod = aBSpl->Period();
        aFlat   = decltype(aFlat)::DownCast(aBSpl->Copy());
        aFlat->SetNotPeriodic();
      }
      const int                         aDeg = aFlat->Degree();
      const int                         aNbP = aFlat->NbPoles();
      const NCollection_Array1<double>& aSeq = aFlat->KnotSequence();
      std::vector<gp_XYZ>               aPoles;
      std::vector<double>               aWeights, aT(aSeq.begin(), aSeq.end());
      for (int i = 1; i <= aNbP; ++i)
      {
        aPoles.push_back(Extrema_CertCC_XYZ(aFlat->Pole(i)));
        aWeights.push_back(aFlat->Weight(i));
      }
      return Extrema_CertCC_SplinePieces(aPoles,
                                         aWeights,
                                         aT,
                                         aDeg,
                                         aPeriod,
                                         theLo,
                                         theHi,
                                         thePieces);
    }
    default:
      return false;
  }
}

// Coaxial circles keep a constant distance where their ranges overlap. The
// global optimizer reports that continuum as parallel; a branch and bound
// cannot discard any cell along it.
inline bool Extrema_CertCC_Coaxial(const gp_Circ& theA, const gp_Circ& theB)
{
  const gp_Ax1& anA = theA.Axis();
  return anA.IsParallel(theB.Axis(), Precision::Angular())
         && gp_Lin(anA).Distance(theB.Location()) <= Precision::Confusion();
}

inline bool Extrema_CertCC_Coaxial(const gp_Circ2d& theA, const gp_Circ2d& theB)
{
  return theA.Location().Distance(theB.Location()) <= Precision::Confusion();
}

struct Extrema_CertCC_Cell
{
  double U0, U1, V0, V1, LB;
  int    A0, A1, B0, B1;

  bool operator<(const Extrema_CertCC_Cell& theOther) const { return LB > theOther.LB; }
};

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePoint>
bool Extrema_CertCC_Perform(const TheCurve1&              theC1,
                            const TheCurve2&              theC2,
                            const double                  theU0,
                            const double                  theU1,
                            const double                  theV0,
                            const double                  theV1,
                            const double                  theTol,
                            const bool                    theIsSingle,
                            NCollection_Sequence<double>& thePoints1,
                            NCollection_Sequence<double>& thePoints2)
{
  using Vec  = typename Extrema_CertCC_VecOf<ThePoint>::Type;
  using Cell = Extrema_CertCC_Cell;
  // Random and test pairs need under 5e4 cells; more means the bounds do not
  // fit the curves, and the global optimizer takes over. A tight cap bounds the
  // work lost on a pair whose minimum is a continuum, which never converges.
  const long aCellCap = 200000;
  const int  aLeafCap = 4096;

  if (TheCurveTool1::GetType(theC1) == GeomAbs_Circle
      && TheCurveTool2::GetType(theC2) == GeomAbs_Circle
      && Extrema_CertCC_Coaxial(TheCurveTool1::Circle(theC1), TheCurveTool2::Circle(theC2)))
  {
    return false;
  }
  std::vector<Extrema_CertCC_Piece> aA, aB;
  if (!Extrema_CertCC_Pieces<TheCurve1, TheCurveTool1>(theC1, theU0, theU1, aA)
      || !Extrema_CertCC_Pieces<TheCurve2, TheCurveTool2>(theC2, theV0, theV1, aB))
  {
    return false;
  }
  const double aTol = std::max(theTol, Precision::PConfusion() * 1.0e-2);

  auto aDist = [&](const double theU, const double theV) {
    return TheCurveTool1::Value(theC1, theU).Distance(TheCurveTool2::Value(theC2, theV));
  };
  double aBestD = RealLast(), aBestU = theU0, aBestV = theV0;
  long   aCells = 0;
  for (const double aU : {theU0, theU1})
  {
    for (const double aV : {theV0, theV1})
    {
      const double aD = aDist(aU, aV);
      if (aD < aBestD)
      {
        aBestD = aD;
        aBestU = aU;
        aBestV = aV;
      }
    }
  }
  // Newton on the squared distance, accepting only steps that lower it. It
  // finds points; only the bounds decide that no lower point exists.
  auto aNewton2d = [&](double& theU, double& theV, double& theD) {
    for (int anIt = 0; anIt < 30; ++anIt)
    {
      ThePoint aPA, aPB;
      Vec      aDA, aDB, anAA, aBB;
      TheCurveTool1::D2(theC1, theU, aPA, aDA, anAA);
      TheCurveTool2::D2(theC2, theV, aPB, aDB, aBB);
      const Vec    aD(aPB, aPA);
      const double aGU = aD.Dot(aDA), aGV = -aD.Dot(aDB);
      const double aHUU = aDA.Dot(aDA) + aD.Dot(anAA), aHVV = aDB.Dot(aDB) - aD.Dot(aBB);
      const double aHUV = -aDA.Dot(aDB);
      const double aDet = aHUU * aHVV - aHUV * aHUV;
      if (!(aHUU > 0.0 && aDet > 0.0))
      {
        return;
      }
      const double aNU = std::min(std::max(theU - (aHVV * aGU - aHUV * aGV) / aDet, theU0), theU1);
      const double aNV = std::min(std::max(theV - (aHUU * aGV - aHUV * aGU) / aDet, theV0), theV1);
      const double aF  = aDist(aNU, aNV);
      if (!(aF < theD))
      {
        return;
      }
      theD = aF;
      theU = aNU;
      theV = aNV;
    }
  };
  // The same along one parameter with the other held at a range end: a
  // minimum on the boundary has a nonzero gradient, which the step above
  // overshoots and clamps short of the end.
  auto aNewton1d = [&](const bool isAlongU, double& theU, double& theV, double& theD) {
    for (int anIt = 0; anIt < 30; ++anIt)
    {
      ThePoint aPA, aPB;
      Vec      aDA, aDB, anAA, aBB;
      TheCurveTool1::D2(theC1, theU, aPA, aDA, anAA);
      TheCurveTool2::D2(theC2, theV, aPB, aDB, aBB);
      const Vec    aD(aPB, aPA);
      const double aG = isAlongU ? aD.Dot(aDA) : -aD.Dot(aDB);
      const double aH = isAlongU ? aDA.Dot(aDA) + aD.Dot(anAA) : aDB.Dot(aDB) - aD.Dot(aBB);
      if (!(aH > 0.0))
      {
        return;
      }
      const double aNU = isAlongU ? std::min(std::max(theU - aG / aH, theU0), theU1) : theU;
      const double aNV = isAlongU ? theV : std::min(std::max(theV - aG / aH, theV0), theV1);
      const double aF  = aDist(aNU, aNV);
      if (!(aF < theD))
      {
        return;
      }
      theD = aF;
      theU = aNU;
      theV = aNV;
    }
  };
  const double aNearU = 1.0e-6 * (theU1 - theU0), aNearV = 1.0e-6 * (theV1 - theV0);
  auto         aNewton = [&](double& theU, double& theV, double& theD) {
    aNewton2d(theU, theV, theD);
    for (const double anEnd : {theU0, theU1})
    {
      if (std::abs(theU - anEnd) <= aNearU)
      {
        double aU = anEnd, aV = theV, aD = aDist(aU, aV);
        aNewton1d(false, aU, aV, aD);
        if (aD <= theD)
        {
          theU = aU;
          theV = aV;
          theD = aD;
        }
      }
    }
    for (const double anEnd : {theV0, theV1})
    {
      if (std::abs(theV - anEnd) <= aNearV)
      {
        double aU = theU, aV = anEnd, aD = aDist(aU, aV);
        aNewton1d(true, aU, aV, aD);
        if (aD <= theD)
        {
          theU = aU;
          theV = aV;
          theD = aD;
        }
      }
    }
  };
  auto aLipschitz =
    [](const std::vector<Extrema_CertCC_Piece>& thePieces, const int theI0, const int theI1) {
      double aL = 0.0;
      for (int i = theI0; i <= theI1; ++i)
      {
        aL = std::max(aL, thePieces[i].L1);
      }
      return aL;
    };
  auto aBound = [&](Cell& theCell) {
    const double aUC = 0.5 * (theCell.U0 + theCell.U1), aVC = 0.5 * (theCell.V0 + theCell.V1);
    const double aHU = 0.5 * (theCell.U1 - theCell.U0), aHV = 0.5 * (theCell.V1 - theCell.V0);
    ThePoint     aPA, aPB;
    Vec          aDA, aDB;
    TheCurveTool1::D1(theC1, aUC, aPA, aDA);
    TheCurveTool2::D1(theC2, aVC, aPB, aDB);
    const Vec    aD(aPB, aPA);
    const double aF = aD.Magnitude();
    ++aCells;
    if (aF < aBestD)
    {
      aBestD = aF;
      aBestU = aUC;
      aBestV = aVC;
    }
    const double aLA = aLipschitz(aA, theCell.A0, theCell.A1);
    const double aLB = aLipschitz(aB, theCell.B0, theCell.B1);
    theCell.LB       = std::max(0.0, aF - aLA * aHU - aLB * aHV);
    if (theCell.A0 == theCell.A1 && theCell.B0 == theCell.B1)
    {
      const double aR   = aF + aLA * aHU + aLB * aHV;
      const double aHuu = 2.0 * (aLA * aLA + aR * aA[theCell.A0].L2);
      const double aHvv = 2.0 * (aLB * aLB + aR * aB[theCell.B0].L2);
      const double aHuv = 2.0 * aLA * aLB;
      const double aG   = aF * aF - std::abs(2.0 * aD.Dot(aDA)) * aHU
                        - std::abs(2.0 * aD.Dot(aDB)) * aHV
                        - 0.5 * (aHuu * aHU * aHU + 2.0 * aHuv * aHU * aHV + aHvv * aHV * aHV);
      if (aG > 0.0)
      {
        theCell.LB = std::max(theCell.LB, std::sqrt(aG));
      }
    }
  };
  auto aSplit = [&](const Cell& theCell, Cell& theL, Cell& theR) {
    const double aWA      = aLipschitz(aA, theCell.A0, theCell.A1) * (theCell.U1 - theCell.U0);
    const double aWB      = aLipschitz(aB, theCell.B0, theCell.B1) * (theCell.V1 - theCell.V0);
    const bool   isMultiA = theCell.A0 < theCell.A1, isMultiB = theCell.B0 < theCell.B1;
    const bool   isSplitU = isMultiA && isMultiB ? aWA >= aWB
                            : isMultiA           ? true
                            : isMultiB           ? false
                                                 : aWA >= aWB;
    theL                  = theCell;
    theR                  = theCell;
    if (isSplitU)
    {
      if (isMultiA)
      {
        const int aM = (theCell.A0 + theCell.A1) / 2;
        theL.A1      = aM;
        theR.A0      = aM + 1;
        theL.U1 = theR.U0 = aA[aM].Hi;
      }
      else
      {
        theL.U1 = theR.U0 = 0.5 * (theCell.U0 + theCell.U1);
      }
    }
    else
    {
      if (isMultiB)
      {
        const int aM = (theCell.B0 + theCell.B1) / 2;
        theL.B1      = aM;
        theR.B0      = aM + 1;
        theL.V1 = theR.V0 = aB[aM].Hi;
      }
      else
      {
        theL.V1 = theR.V0 = 0.5 * (theCell.V0 + theCell.V1);
      }
    }
  };

  const Cell aRootTemplate{theU0,
                           theU1,
                           theV0,
                           theV1,
                           0.0,
                           0,
                           static_cast<int>(aA.size()) - 1,
                           0,
                           static_cast<int>(aB.size()) - 1};
  Cell       aRoot = aRootTemplate;
  aBound(aRoot);
  std::priority_queue<Cell> anOpen;
  anOpen.push(aRoot);
  long aNextNewton = 1;
  while (!anOpen.empty() && anOpen.top().LB < aBestD - aTol)
  {
    if (aCells > aCellCap)
    {
      return false;
    }
    if (aCells >= aNextNewton)
    {
      aNewton(aBestU, aBestV, aBestD);
      aNextNewton *= 4;
    }
    const Cell aCell = anOpen.top();
    anOpen.pop();
    Cell aL, aR;
    aSplit(aCell, aL, aR);
    aBound(aL);
    aBound(aR);
    if (aL.LB < aBestD - aTol)
    {
      anOpen.push(aL);
    }
    if (aR.LB < aBestD - aTol)
    {
      anOpen.push(aR);
    }
  }
  aNewton(aBestU, aBestV, aBestD);

  thePoints1.Clear();
  thePoints2.Clear();
  if (theIsSingle)
  {
    thePoints1.Append(aBestU);
    thePoints2.Append(aBestV);
    return true;
  }

  // Every point within theTol of the minimum. Cells that may hold one are cut
  // down to leaves, Newton runs from each leaf, and touching leaves are one
  // solution.
  const double      aThreshold = aBestD + aTol;
  const double      aLeafWidth = std::max(1.0e-3
                                       * (aLipschitz(aA, 0, aRoot.A1) * (theU1 - theU0)
                                          + aLipschitz(aB, 0, aRoot.B1) * (theV1 - theV0)),
                                     4.0 * aTol);
  std::vector<Cell> aStack(1, aRoot), aLeaves;
  while (!aStack.empty())
  {
    if (aCells > aCellCap || static_cast<int>(aLeaves.size()) > aLeafCap)
    {
      return false;
    }
    const Cell aCell = aStack.back();
    aStack.pop_back();
    if (aCell.LB > aThreshold)
    {
      continue;
    }
    const bool isLeaf =
      aCell.A0 == aCell.A1 && aCell.B0 == aCell.B1
      && aA[aCell.A0].L1 * (aCell.U1 - aCell.U0) + aB[aCell.B0].L1 * (aCell.V1 - aCell.V0)
           <= aLeafWidth;
    if (isLeaf)
    {
      aLeaves.push_back(aCell);
      continue;
    }
    Cell aL, aR;
    aSplit(aCell, aL, aR);
    aBound(aL);
    aBound(aR);
    aStack.push_back(aL);
    aStack.push_back(aR);
  }

  const int        aNbLeaves = static_cast<int>(aLeaves.size());
  std::vector<int> aParent(aNbLeaves);
  for (int i = 0; i < aNbLeaves; ++i)
  {
    aParent[i] = i;
  }
  std::function<int(int)> aFind = [&](const int theI) {
    return aParent[theI] == theI ? theI : (aParent[theI] = aFind(aParent[theI]));
  };
  for (int i = 0; i < aNbLeaves; ++i)
  {
    for (int j = i + 1; j < aNbLeaves; ++j)
    {
      const Cell& aP = aLeaves[i];
      const Cell& aQ = aLeaves[j];
      if (aP.U0 <= aQ.U1 && aQ.U0 <= aP.U1 && aP.V0 <= aQ.V1 && aQ.V0 <= aP.V1)
      {
        aParent[aFind(i)] = aFind(j);
      }
    }
  }

  struct Point
  {
    double U, V, D;
  };

  struct Group
  {
    Point Best{0.0, 0.0, RealLast()}, LoU, HiU, LoV, HiV;
    bool  IsAny = false;
  };

  std::vector<Group> aGroups(aNbLeaves);
  for (int i = 0; i < aNbLeaves; ++i)
  {
    const Cell& aCell = aLeaves[i];
    double      aU = 0.5 * (aCell.U0 + aCell.U1), aV = 0.5 * (aCell.V0 + aCell.V1);
    double      aD = aDist(aU, aV);
    if (aCell.U0 <= aBestU && aBestU <= aCell.U1 && aCell.V0 <= aBestV && aBestV <= aCell.V1)
    {
      aU = aBestU;
      aV = aBestV;
      aD = aBestD;
    }
    aNewton(aU, aV, aD);
    Group&      aG = aGroups[aFind(i)];
    const Point aP{aU, aV, aD};
    if (aD < aG.Best.D)
    {
      aG.Best = aP;
    }
    if (aD > aThreshold)
    {
      continue;
    }
    if (!aG.IsAny)
    {
      aG = Group{aP, aP, aP, aP, aP, true};
      continue;
    }
    aG.LoU = aU < aG.LoU.U ? aP : aG.LoU;
    aG.HiU = aU > aG.HiU.U ? aP : aG.HiU;
    aG.LoV = aV < aG.LoV.V ? aP : aG.LoV;
    aG.HiV = aV > aG.HiV.V ? aP : aG.HiV;
  }
  // A stretch of points at the minimum, such as two curves running together,
  // is reported by its two ends, as the global optimizer reports it.
  const double aLA = aLipschitz(aA, 0, aRoot.A1), aLB = aLipschitz(aB, 0, aRoot.B1);
  auto         aSpan = [&](const Point& theP, const Point& theQ) {
    return aLA * std::abs(theP.U - theQ.U) + aLB * std::abs(theP.V - theQ.V);
  };
  // Moves an end outward along the stretch while the distance stays within a
  // small part of theTol of the minimum, so it does not stop at the centre of
  // its leaf and does not report a distance worse than the minimum found.
  const double anEndThreshold = aBestD + 1.0e-3 * aTol;
  auto         aPushOut       = [&](const Point& theEnd, const Point& theOther) {
    const double aLen = aSpan(theEnd, theOther);
    const double aDU = (theEnd.U - theOther.U) / aLen, aDV = (theEnd.V - theOther.V) / aLen;
    auto         anAt = [&](const double theS) {
      const double aU = std::min(std::max(theEnd.U + aDU * theS, theU0), theU1);
      const double aV = std::min(std::max(theEnd.V + aDV * theS, theV0), theV1);
      return Point{aU, aV, aDist(aU, aV)};
    };
    const Point aFar = anAt(2.0 * aLeafWidth);
    if (aFar.D <= anEndThreshold)
    {
      return aFar;
    }
    double aIn = 0.0, anOut = 2.0 * aLeafWidth;
    Point  aBest = theEnd;
    for (int anIt = 0; anIt < 60 && anOut - aIn > 1.0e-6 * aTol; ++anIt)
    {
      const double aMid = 0.5 * (aIn + anOut);
      const Point  aP   = anAt(aMid);
      if (aP.D <= anEndThreshold)
      {
        aIn   = aMid;
        aBest = aP;
      }
      else
      {
        anOut = aMid;
      }
    }
    return aBest;
  };
  for (int i = 0; i < aNbLeaves; ++i)
  {
    const Group& aG = aGroups[i];
    if (aFind(i) != i || !aG.IsAny)
    {
      continue;
    }
    const bool   isAlongU = aSpan(aG.LoU, aG.HiU) >= aSpan(aG.LoV, aG.HiV);
    const Point& aLo      = isAlongU ? aG.LoU : aG.LoV;
    const Point& aHi      = isAlongU ? aG.HiU : aG.HiV;
    if (aSpan(aLo, aHi) > 2.0 * aLeafWidth)
    {
      for (const Point& anEnd : {aPushOut(aLo, aHi), aPushOut(aHi, aLo)})
      {
        thePoints1.Append(anEnd.U);
        thePoints2.Append(anEnd.V);
      }
    }
    else
    {
      thePoints1.Append(aG.Best.U);
      thePoints2.Append(aG.Best.V);
    }
  }
  if (thePoints1.IsEmpty())
  {
    thePoints1.Append(aBestU);
    thePoints2.Append(aBestV);
  }
  return true;
}
} // namespace

//=================================================================================================
// Implementation
//=================================================================================================

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
Extrema_GGenExtCC<TheCurve1, TheCurveTool1, TheCurve2, TheCurveTool2, ThePOnC, ThePoint, TheExtPC>::
  Extrema_GGenExtCC()
    : myIsFindSingleSolution(false),
      myParallel(false),
      myCurveMinTol(Precision::PConfusion()),
      myLowBorder(1, 2),
      myUppBorder(1, 2),
      myDone(false)
{
  myC[0] = myC[1] = nullptr;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
Extrema_GGenExtCC<TheCurve1, TheCurveTool1, TheCurve2, TheCurveTool2, ThePOnC, ThePoint, TheExtPC>::
  Extrema_GGenExtCC(const TheCurve1& theC1, const TheCurve2& theC2)
    : myIsFindSingleSolution(false),
      myParallel(false),
      myCurveMinTol(Precision::PConfusion()),
      myLowBorder(1, 2),
      myUppBorder(1, 2),
      myDone(false)
{
  myC[0]         = (void*)&theC1;
  myC[1]         = (void*)&theC2;
  myLowBorder(1) = theC1.FirstParameter();
  myLowBorder(2) = theC2.FirstParameter();
  myUppBorder(1) = theC1.LastParameter();
  myUppBorder(2) = theC2.LastParameter();
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
Extrema_GGenExtCC<TheCurve1, TheCurveTool1, TheCurve2, TheCurveTool2, ThePOnC, ThePoint, TheExtPC>::
  Extrema_GGenExtCC(const TheCurve1& theC1,
                    const TheCurve2& theC2,
                    const double     theUinf,
                    const double     theUsup,
                    const double     theVinf,
                    const double     theVsup)
    : myIsFindSingleSolution(false),
      myParallel(false),
      myCurveMinTol(Precision::PConfusion()),
      myLowBorder(1, 2),
      myUppBorder(1, 2),
      myDone(false)
{
  myC[0]         = (void*)&theC1;
  myC[1]         = (void*)&theC2;
  myLowBorder(1) = theUinf;
  myLowBorder(2) = theVinf;
  myUppBorder(1) = theUsup;
  myUppBorder(2) = theVsup;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
void Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::SetParams(const TheCurve1& theC1,
                                            const TheCurve2& theC2,
                                            const double     theUinf,
                                            const double     theUsup,
                                            const double     theVinf,
                                            const double     theVsup)
{
  myC[0]         = (void*)&theC1;
  myC[1]         = (void*)&theC2;
  myLowBorder(1) = theUinf;
  myLowBorder(2) = theVinf;
  myUppBorder(1) = theUsup;
  myUppBorder(2) = theVsup;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
void Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::SetTolerance(const double theTol)
{
  myCurveMinTol = theTol;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
void Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::SetSingleSolutionFlag(const bool theFlag)
{
  myIsFindSingleSolution = theFlag;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
bool Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::GetSingleSolutionFlag() const
{
  return myIsFindSingleSolution;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
void Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::Perform()
{
  myDone     = false;
  myParallel = false;

  TheCurve1& C1 = *(TheCurve1*)myC[0];
  TheCurve2& C2 = *(TheCurve2*)myC[1];

  if (Extrema_CertCC_Perform<TheCurve1, TheCurveTool1, TheCurve2, TheCurveTool2, ThePoint>(
        C1,
        C2,
        myLowBorder(1),
        myUppBorder(1),
        myLowBorder(2),
        myUppBorder(2),
        myCurveMinTol,
        myIsFindSingleSolution,
        myPoints1,
        myPoints2))
  {
    myDone = true;
    return;
  }

  int           aNbInter[2];
  GeomAbs_Shape aContinuity = GeomAbs_C2;
  aNbInter[0]               = C1.NbIntervals(aContinuity);
  aNbInter[1]               = C2.NbIntervals(aContinuity);

  if (aNbInter[0] * aNbInter[1] > 100)
  {
    aContinuity = GeomAbs_C1;
    aNbInter[0] = C1.NbIntervals(aContinuity);
    aNbInter[1] = C2.NbIntervals(aContinuity);
  }

  double       anL[2];
  int          indmax = -1, indmin = -1;
  const double mult = 20.;
  if (!(Precision::IsInfinite(C1.FirstParameter()) || Precision::IsInfinite(C1.LastParameter())
        || Precision::IsInfinite(C2.FirstParameter()) || Precision::IsInfinite(C2.LastParameter())))
  {
    anL[0] = GCPnts_AbscissaPoint::Length(C1);
    anL[1] = GCPnts_AbscissaPoint::Length(C2);
    if (anL[0] / aNbInter[0] > mult * anL[1] / aNbInter[1])
    {
      indmax = 0;
      indmin = 1;
    }
    else if (anL[1] / aNbInter[1] > mult * anL[0] / aNbInter[0])
    {
      indmax = 1;
      indmin = 0;
    }
  }
  int aNbIntOpt = 0;
  if (indmax >= 0)
  {
    aNbIntOpt = RealToInt(anL[indmax] * aNbInter[indmin] / anL[indmin] / (mult / 4.)) + 1;
    if (aNbIntOpt > 100 || aNbIntOpt < aNbInter[indmax])
    {
      indmax = -1;
    }
    else
    {
      if (aNbIntOpt * aNbInter[indmin] > 100)
      {
        aNbIntOpt = 100 / aNbInter[indmin];
        if (aNbIntOpt < aNbInter[indmax])
        {
          indmax = -1;
        }
      }
    }
  }

  occ::handle<NCollection_HArray1<double>> anIntervals1 =
    new NCollection_HArray1<double>(1, aNbInter[0] + 1);
  occ::handle<NCollection_HArray1<double>> anIntervals2 =
    new NCollection_HArray1<double>(1, aNbInter[1] + 1);
  C1.Intervals(anIntervals1->ChangeArray1(), aContinuity);
  C2.Intervals(anIntervals2->ChangeArray1(), aContinuity);
  if (indmax >= 0)
  {
    if (indmax == 0)
    {
      Extrema_GGenExtCC_ChangeIntervals(anIntervals1, aNbIntOpt);
      aNbInter[0] = anIntervals1->Length() - 1;
    }
    else
    {
      Extrema_GGenExtCC_ChangeIntervals(anIntervals2, aNbIntOpt);
      aNbInter[1] = anIntervals2->Length() - 1;
    }
  }
  if (C1.IsClosed() && aNbInter[0] == 1)
  {
    Extrema_GGenExtCC_ChangeIntervals(anIntervals1, 3);
    aNbInter[0] = anIntervals1->Length() - 1;
  }
  if (C2.IsClosed() && aNbInter[1] == 1)
  {
    Extrema_GGenExtCC_ChangeIntervals(anIntervals2, 3);
    aNbInter[1] = anIntervals2->Length() - 1;
  }

  const double aMaxLC   = 10000.;
  double       aLC      = 100.0;
  const double aMaxDer1 = 1.0 / C1.Resolution(1.0);
  const double aMaxDer2 = 1.0 / C2.Resolution(1.0);
  double       aMaxDer  = std::max(aMaxDer1, aMaxDer2) * M_SQRT2;
  if (aLC > aMaxDer)
    aLC = aMaxDer;

  bool         isConstLockedFlag = false;
  const double aCR               = 0.001;
  if (aMaxDer1 / aMaxDer < aCR || aMaxDer2 / aMaxDer < aCR)
  {
    isConstLockedFlag = true;
  }
  if (aMaxDer > aMaxLC)
  {
    aLC               = aMaxLC;
    isConstLockedFlag = true;
  }
  if (C1.GetType() == GeomAbs_Line)
  {
    aMaxDer = 1.0 / C2.Resolution(1.0);
    if (aLC > aMaxDer)
    {
      isConstLockedFlag = true;
      aLC               = aMaxDer;
    }
  }
  if (C2.GetType() == GeomAbs_Line)
  {
    aMaxDer = 1.0 / C1.Resolution(1.0);
    if (aLC > aMaxDer)
    {
      isConstLockedFlag = true;
      aLC               = aMaxDer;
    }
  }

  Extrema_GlobOptFuncCCC2 aFunc(C1, C2);
  if (aLC < aMaxLC || aMaxDer > aMaxLC)
  {
    math_Vector aT(1, 2), aG(1, 2);
    double      aF, aMaxG = 0.;
    double      t1, t2, dt1, dt2;
    int         n1 = 21, n2 = 21, i1, i2;
    dt1 = (C1.LastParameter() - C1.FirstParameter()) / (n1 - 1);
    dt2 = (C2.LastParameter() - C2.FirstParameter()) / (n2 - 1);
    for (i1 = 1, t1 = C1.FirstParameter(); i1 <= n1; ++i1, t1 += dt1)
    {
      aT(1) = t1;
      for (i2 = 1, t2 = C2.FirstParameter(); i2 <= n2; ++i2, t2 += dt2)
      {
        aT(2) = t2;
        aFunc.Values(aT, aF, aG);
        double aMod = aG(1) * aG(1) + aG(2) * aG(2);
        aMaxG       = std::max(aMaxG, aMod);
      }
    }
    aMaxG = std::sqrt(aMaxG);
    if (aMaxG > aMaxDer)
    {
      aLC               = std::min(aMaxG, aMaxLC);
      isConstLockedFlag = true;
    }
    if (aMaxG > 100. * aMaxLC)
    {
      aLC               = 100. * aMaxLC;
      isConstLockedFlag = true;
    }
    else if (aMaxG < 0.1 * aMaxDer)
    {
      isConstLockedFlag = true;
    }
  }
  math_GlobOptMin aFinder(&aFunc, myLowBorder, myUppBorder, aLC);
  aFinder.SetLipConstState(isConstLockedFlag);
  aFinder.SetContinuity(aContinuity == GeomAbs_C2 ? 2 : 1);
  double aDiscTol  = 1.0e-2;
  double aValueTol = 1.0e-2;
  double aSameTol  = myCurveMinTol / (aDiscTol);
  aFinder.SetTol(aDiscTol, aSameTol);
  aFinder.SetFunctionalMinimalValue(0.0);

  const double aCellSize = std::max(std::max(anIntervals1->Last() - anIntervals1->First(),
                                             anIntervals2->Last() - anIntervals2->First())
                                      * Precision::PConfusion() / (2.0 * M_SQRT2),
                                    Precision::PConfusion());
  Extrema_GGenExtCC_PointsInspector                         anInspector(aCellSize);
  NCollection_CellFilter<Extrema_GGenExtCC_PointsInspector> aFilter(aCellSize);
  NCollection_DynamicArray<gp_XY>                           aPnts;

  int         i, j, k;
  math_Vector aFirstBorderInterval(1, 2);
  math_Vector aSecondBorderInterval(1, 2);
  double      aF     = RealLast();
  double      aCurrF = RealLast();
  for (i = 1; i <= aNbInter[0]; i++)
  {
    for (j = 1; j <= aNbInter[1]; j++)
    {
      aFirstBorderInterval(1)  = anIntervals1->Value(i);
      aFirstBorderInterval(2)  = anIntervals2->Value(j);
      aSecondBorderInterval(1) = anIntervals1->Value(i + 1);
      aSecondBorderInterval(2) = anIntervals2->Value(j + 1);

      aFinder.SetLocalParams(aFirstBorderInterval, aSecondBorderInterval);
      aFinder.Perform(GetSingleSolutionFlag());

      aCurrF = aFinder.GetF();
      if (aCurrF >= aF + aSameTol * aValueTol)
      {
        continue;
      }

      if (aCurrF > aF - aSameTol * aValueTol)
      {
        if (aCurrF < aF)
          aF = aCurrF;
      }
      else
      {
        aF = aCurrF;
        aFilter.Reset(aCellSize);
        aPnts.Clear();
      }

      math_Vector sol(1, 2);
      for (k = 1; k <= aFinder.NbExtrema(); k++)
      {
        aFinder.Points(k, sol);
        gp_XY aPnt2d(sol(1), sol(2));

        gp_XY aXYmin = anInspector.Shift(aPnt2d, -aCellSize);
        gp_XY aXYmax = anInspector.Shift(aPnt2d, aCellSize);

        anInspector.ClearFind();
        anInspector.SetCurrent(aPnt2d);
        aFilter.Inspect(aXYmin, aXYmax, anInspector);
        if (!anInspector.isFind())
        {
          aFilter.Add(aPnt2d, aPnt2d);
          aPnts.Append(gp_XY(sol(1), sol(2)));
        }
      }
    }
  }

  const int aNbSol = aPnts.Length();
  if (aNbSol == 0)
  {
    myDone = false;
    return;
  }

  myDone = true;

  if (aNbSol == 1)
  {
    const gp_XY& aSol = aPnts.First();
    myPoints1.Append(aSol.X());
    myPoints2.Append(aSol.Y());
    return;
  }

  std::sort(aPnts.begin(), aPnts.end(), Extrema_GGenExtCC_comp);

  NCollection_List<int> aSolutions;

  bool bSaveSolution       = true;
  bool bDirsCoinside       = true;
  bool bDifferentSolutions = false;

  bool        isParallel = true;
  double      aVal       = 0.0;
  math_Vector aVec(1, 2, 0.0);

  for (int anIdx = 0; anIdx < aNbSol - 1; anIdx++)
  {
    const gp_XY& aCurrent = aPnts(anIdx);
    const gp_XY& aNext    = aPnts(anIdx + 1);

    aVec(1) = (aCurrent.X() + aNext.X()) * 0.5;
    aVec(2) = (aCurrent.Y() + aNext.Y()) * 0.5;

    aFunc.Value(aVec, aVal);
    if (std::abs(aVal - aF) < Precision::Confusion())
    {
      if (bSaveSolution)
      {
        aSolutions.Append(anIdx);
        bSaveSolution = false;
      }
    }
    else
    {
      isParallel = false;
      aSolutions.Append(anIdx);
      bSaveSolution = true;
    }

    if (!bDifferentSolutions)
    {
      if (aNext.X() > aCurrent.X())
      {
        if (aNext.Y() > aCurrent.Y())
        {
          bDifferentSolutions = true;
          bDirsCoinside       = true;
        }
        else if (aNext.Y() < aCurrent.Y())
        {
          bDifferentSolutions = true;
          bDirsCoinside       = false;
        }
      }
    }
  }
  aSolutions.Append(aNbSol - 1);

  if (!bDifferentSolutions)
    isParallel = false;

  if (isParallel)
  {
    double aT1[2] = {myLowBorder(1), myUppBorder(1)};
    double aT2[2] = {bDirsCoinside ? myLowBorder(2) : myUppBorder(2),
                     bDirsCoinside ? myUppBorder(2) : myLowBorder(2)};

    TheExtPC anExtPC1, anExtPC2;
    anExtPC1.Initialize(C1, myLowBorder(1), myUppBorder(1));
    anExtPC2.Initialize(C2, myLowBorder(2), myUppBorder(2));

    for (int iT = 0; isParallel && (iT < 2); ++iT)
    {
      double aDist1 =
        Extrema_GGenExtCC_ProjPOnC<TheCurve2, TheExtPC, ThePoint>(C1.Value(aT1[iT]), anExtPC2);
      double aDist2 =
        Extrema_GGenExtCC_ProjPOnC<TheCurve1, TheExtPC, ThePoint>(C2.Value(aT2[iT]), anExtPC1);
      isParallel = (std::abs(std::min(aDist1, aDist2) - aF * aF) < Precision::Confusion());
    }
  }

  if (isParallel)
  {
    const gp_XY& aSol = aPnts.First();
    myPoints1.Append(aSol.X());
    myPoints2.Append(aSol.Y());
    myParallel = true;
  }
  else
  {
    NCollection_List<int>::Iterator aItSol(aSolutions);
    for (; aItSol.More(); aItSol.Next())
    {
      const gp_XY& aSol = aPnts(aItSol.Value());
      myPoints1.Append(aSol.X());
      myPoints2.Append(aSol.Y());
    }
  }
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
bool Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::IsDone() const
{
  return myDone;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
bool Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::IsParallel() const
{
  if (!IsDone())
    throw StdFail_NotDone();
  return myParallel;
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
int Extrema_GGenExtCC<TheCurve1,
                      TheCurveTool1,
                      TheCurve2,
                      TheCurveTool2,
                      ThePOnC,
                      ThePoint,
                      TheExtPC>::NbExt() const
{
  if (!IsDone())
    throw StdFail_NotDone();
  return myPoints1.Length();
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
double Extrema_GGenExtCC<TheCurve1,
                         TheCurveTool1,
                         TheCurve2,
                         TheCurveTool2,
                         ThePOnC,
                         ThePoint,
                         TheExtPC>::SquareDistance(const int theN) const
{
  if (theN < 1 || theN > NbExt())
  {
    throw Standard_OutOfRange();
  }

  return TheCurveTool1::Value(*((TheCurve1*)myC[0]), myPoints1(theN))
    .SquareDistance(TheCurveTool2::Value(*((TheCurve2*)myC[1]), myPoints2(theN)));
}

//=================================================================================================

template <typename TheCurve1,
          typename TheCurveTool1,
          typename TheCurve2,
          typename TheCurveTool2,
          typename ThePOnC,
          typename ThePoint,
          typename TheExtPC>
void Extrema_GGenExtCC<TheCurve1,
                       TheCurveTool1,
                       TheCurve2,
                       TheCurveTool2,
                       ThePOnC,
                       ThePoint,
                       TheExtPC>::Points(const int theN, ThePOnC& theP1, ThePOnC& theP2) const
{
  if (theN < 1 || theN > NbExt())
  {
    throw Standard_OutOfRange();
  }

  theP1.SetValues(myPoints1(theN), TheCurveTool1::Value(*((TheCurve1*)myC[0]), myPoints1(theN)));
  theP2.SetValues(myPoints2(theN), TheCurveTool2::Value(*((TheCurve2*)myC[1]), myPoints2(theN)));
}

#endif // _Extrema_GGenExtCC_HeaderFile

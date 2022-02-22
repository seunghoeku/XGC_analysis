//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtkm_exec_celllocatorxgcgrid_h
#define vtkm_exec_celllocatorxgcgrid_h

#include <vtkm/Bounds.h>
#include <vtkm/Math.h>
#include <vtkm/TopologyElementTag.h>
#include <vtkm/Types.h>
#include <vtkm/VecFromPortalPermute.h>

#include <vtkm/cont/ArrayHandleXGCCoordinates.h>
#include <vtkm/cont/CellSetExtrude.h>
#include <vtkm/cont/CellSetSingleType.h>

#include <vtkm/exec/CellInside.h>
#include <vtkm/exec/CellLocatorTwoLevel.h>
#include <vtkm/exec/ConnectivityExtrude.h>
#include <vtkm/exec/ParametricCoordinates.h>

namespace vtkm
{

namespace exec
{

using FloatVec3 = vtkm::Vec3f;
using CellLocatorType =
  vtkm::cont::CellSetSingleType<>::ExecConnectivityType<vtkm::TopologyElementTagCell,
                                                        vtkm::TopologyElementTagPoint>;

using CoordsPortalType = vtkm::internal::ArrayPortalXGCCoordinates<
  vtkm::internal::ArrayPortalBasicRead<vtkm::FloatDefault>>;

class VTKM_ALWAYS_EXPORT CellLocatorXGCGrid
{
public:
  VTKM_CONT
  CellLocatorXGCGrid(const vtkm::exec::ConnectivityExtrude& conn,
                     const CoordsPortalType& coords,
                     const vtkm::exec::CellLocatorTwoLevel<CellLocatorType>& planeLocator,
                     const vtkm::Id& numPlanes,
                     const vtkm::Id& cellsPerPlane,
                     const bool& useCylindrical)
    : CellsPerPlane(cellsPerPlane)
    , Connectivity(conn)
    , Coords(coords)
    , NumPlanes(numPlanes)
    , PlaneLocator(planeLocator)
    , ThetaSpacing(vtkm::TwoPi() / numPlanes)
    , UseCylindrical(useCylindrical)
  {
  }

  VTKM_EXEC
  vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                           vtkm::Id& cellId,
                           vtkm::Vec3f& parametric) const
  {
    vtkm::Vec3f cylPt;
    if (this->UseCylindrical)
    {
      cylPt = point;
    }
    else
    {
      vtkm::FloatDefault x = point[0];
      vtkm::FloatDefault y = point[1];
      vtkm::FloatDefault z = point[2];

      vtkm::FloatDefault r = vtkm::Sqrt(x * x + y * y);
      vtkm::FloatDefault theta = vtkm::ATan2(y, x);
      if (theta < 0)
        theta += vtkm::TwoPi();
      cylPt[0] = r;
      cylPt[1] = theta;
      cylPt[2] = z;
    }
    vtkm::FloatDefault theta = cylPt[1];

    vtkm::Vec3f cylPt2D(cylPt[0], cylPt[2], 0);
    //std::cout << "FindCell: " << point << " --> " << cylPt << " 2d= "<<cylPt2D<<" theta= "<<(cylPt[1]*57.2958)<<std::endl;

    vtkm::Id cid = -1;
    auto res = this->PlaneLocator.FindCell(cylPt2D, cid, parametric);
    //std::cout<<"     plane cid= "<<cid<<std::endl;

    if (res != vtkm::ErrorCode::Success)
      return res;

    vtkm::Id planeIdx = vtkm::Floor(theta / this->ThetaSpacing);

    if (planeIdx > 0)
      cid += (planeIdx * this->CellsPerPlane);

    auto indices = this->Connectivity.GetIndices(cid);
    //std::cout<<"CID= "<<cid<<" planeIdx= "<<planeIdx<<std::endl;

    auto pts = vtkm::make_VecFromPortalPermute(&indices, this->Coords);
    //for (int i = 0; i < 6; i++)
      //std::cout << "Pt_" << i << " idx= " << indices[i] << " pt= "<<pts[i]<<std::endl;

    //do the wedge tests in R,Theta,Z space.
    /*
    cylPt[0] = r;
    cylPt[1] = theta;
    cylPt[2] = z;
    */

    vtkm::Vec<vtkm::Vec3f, 6> cylVerts, cylVertsDeg;
    for (int i = 0; i < 6; i++)
    {
      vtkm::Vec3f p = pts[i];
      if (this->UseCylindrical)
      {
        cylVerts[i] = p;
      }
      else
      {
        vtkm::FloatDefault r = vtkm::Sqrt(p[0]*p[0] + p[1]*p[1]);
        vtkm::FloatDefault pTheta = vtkm::ATan2(p[1], p[0]);
        if (pTheta < 0) pTheta += vtkm::TwoPi();
        cylVerts[i][0] = r;
        cylVerts[i][1] = pTheta;
        cylVerts[i][2] = p[2];
      }
//      if (cylVerts[i][1] < 0)
//        cylVerts[i][1] += vtkm::TwoPi();
      cylVertsDeg[i] = cylVerts[i];
      cylVertsDeg[i][1] = p[1] * 57.2958;
    }

    //Wrap around. Last plane at 0, so do a wraparound.
    //Need to handle the other case: last plane at 2pi
    if (this->UseCylindrical && planeIdx == (this->NumPlanes-1))
    {
      //Last plane should be at 0,
      VTKM_ASSERT(cylVerts[0][1] > cylVerts[3][1]);
      VTKM_ASSERT(cylVerts[1][1] > cylVerts[4][1]);
      VTKM_ASSERT(cylVerts[2][1] > cylVerts[5][1]);
      cylVerts[3][1] += vtkm::TwoPi();
      cylVerts[4][1] += vtkm::TwoPi();
      cylVerts[5][1] += vtkm::TwoPi();
    }

    auto xx = cylPt;
    xx[1] *=  57.2958;
    //std::cout<<"****** cylPt= "<<xx<<std::endl;
    //for (int i = 0; i < 6; i++)
    //std::cout << "CPt_" << i << " idx= " << indices[i] << " pt= "<<cylVertsDeg[i]<<std::endl;


    FloatVec3 pc;
    bool inside;
    //VTKM_RETURN_ON_ERROR(this->PointInsideCell(point, pts, pc, inside));
    VTKM_RETURN_ON_ERROR(this->PointInsideCell(cylPt, cylVerts, pc, inside));
    if (inside)
    {
      cellId = cid;
      parametric = pc;

      //std::cout << " *** cellId= " << cellId << " p: " << parametric << std::endl;
      //std::cout << "   planeIdx= " << planeIdx << " cellId= " << cellId << std::endl;
      //std::cout << "   parametric= " << parametric << std::endl;
      //std::cout<<std::endl<<std::endl;
      //std::cout<<"************************ FOUND!!!!!!!!!!!"<<std::endl;

      return vtkm::ErrorCode::Success;
    }

    //std::cout<<"************************ NOT FOUND"<<std::endl;
    printf("XGC CellNotFound %d\n", __LINE__);
    return vtkm::ErrorCode::CellNotFound;
  }

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorXGCGrid* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorXGCGrid* operator->() const { return this; }

private:
  template <typename PointsVecType>
  VTKM_EXEC vtkm::Bounds ComputeCellBounds(const PointsVecType& points) const
  {
    using CoordsType = typename vtkm::VecTraits<PointsVecType>::ComponentType;
    CoordsType minp = points[0], maxp = points[0];

    for (vtkm::IdComponent i = 1; i < 6; i++)
    {
      minp = vtkm::Min(minp, points[i]);
      maxp = vtkm::Max(maxp, points[i]);
    }

    return { FloatVec3(minp), FloatVec3(maxp) };
  }

  VTKM_EXEC bool InBounds(const FloatVec3& point,
                          const vtkm::Bounds& bounds,
                          const vtkm::FloatDefault& eps) const
  {
#define isBetween(A, B, C) ( ((A-B) > -eps) && ((A-C) < eps) )

    //std::cout<<"   InBounds: "<<point<<" "<<bounds;
    if (isBetween(point[0], bounds.X.Min, bounds.X.Max) &&
        isBetween(point[1], bounds.Y.Min, bounds.Y.Max) &&
        isBetween(point[2], bounds.Z.Min, bounds.Z.Max))
    {
      //std::cout<<" -->INSIDE"<<std::endl;
      return true;
    }
    //std::cout<<" -->OUTSIDE"<<std::endl;
    return false;
  }


  template <typename CoordsType>
  VTKM_EXEC vtkm::ErrorCode PointInsideCell(FloatVec3 point,
                                            CoordsType cellPoints,
                                            FloatVec3& parametricCoordinates,
                                            bool& inside) const
  {
    vtkm::Bounds bounds = this->ComputeCellBounds(cellPoints);

    inside = false;
    vtkm::FloatDefault eps = 1e-6;
    if (this->InBounds(point, bounds, eps))
    {
      VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
        cellPoints, point, vtkm::CellShapeTagWedge{}, parametricCoordinates));
      //std::cout<<"   PARAMETRIC: "<<parametricCoordinates<<std::endl;
      inside = vtkm::exec::CellInside(parametricCoordinates, vtkm::CellShapeTagWedge{});
    }

    // Return success error code even point is not inside this cell
    return vtkm::ErrorCode::Success;
  }

  template <typename CoordsType>
  VTKM_EXEC vtkm::ErrorCode PointInsideCellCyl(FloatVec3 point,
                                               CoordsType cellPoints,
                                               FloatVec3& parametricCoordinates,
                                               bool& inside) const
  {
    vtkm::Bounds bounds;// = this->ComputeCellBounds(cellPoints);

    inside = false;
    vtkm::FloatDefault eps = 1e-6;
    if (this->InBounds(point, bounds, eps))
    {
      VTKM_RETURN_ON_ERROR(vtkm::exec::WorldCoordinatesToParametricCoordinates(
        cellPoints, point, vtkm::CellShapeTagWedge{}, parametricCoordinates));
      //std::cout<<"   PARAMETRIC: "<<parametricCoordinates<<std::endl;
      inside = vtkm::exec::CellInside(parametricCoordinates, vtkm::CellShapeTagWedge{});
    }

    // Return success error code even point is not inside this cell
    return vtkm::ErrorCode::Success;
  }

/*
  VTKM_EXEC
  vtkm::ErrorCode FindCellCylindrical(const vtkm::Vec3f& point,
                                      vtkm::Id& cellId,
                                      vtkm::Vec3f& parametric) const
  {
  }

  VTKM_EXEC
  vtkm::ErrorCode FindCellCartesian(const vtkm::Vec3f& point,
                                    vtkm::Id& cellId,
                                    vtkm::Vec3f& parametric) const
  {
  }
*/

  vtkm::Id CellsPerPlane;
  vtkm::exec::ConnectivityExtrude Connectivity;
  CoordsPortalType Coords;
  vtkm::Id NumPlanes;
  vtkm::exec::CellLocatorTwoLevel<CellLocatorType> PlaneLocator;
  vtkm::FloatDefault ThetaSpacing;
  bool UseCylindrical;
};
}
}

#endif //vtkm_exec_celllocatorxgcgrid_h

//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_exec_CellLocatorMultiplexer_h
#define vtk_m_exec_CellLocatorMultiplexer_h

#include <vtkm/ErrorCode.h>
#include <vtkm/TypeList.h>

#include <vtkm/exec/internal/Variant.h>

namespace vtkm
{
namespace exec
{

using DimensionType = vtkm::Int16;
using DimVec3 = vtkm::Vec<DimensionType, 3>;
using FloatVec3 = vtkm::Vec3f;
struct Grid
{
  DimVec3 Dimensions;
  // Bug in CUDA 9.2 where having this gap for alignment was for some reason setting garbage
  // in a union with other cell locators (or perhaps not properly copying data). This appears
  // to be fixed by CUDA 10.2.
  DimensionType Padding;
  FloatVec3 Origin;
  FloatVec3 BinSize;
};

namespace detail
{
template <typename T>
using ReadPortal = typename vtkm::cont::ArrayHandle<T>::ReadPortalType;

struct FindCellFunctor
{
  template <typename Locator>
  VTKM_EXEC vtkm::ErrorCode operator()(Locator&& locator,
                                       const vtkm::Vec3f& point,
                                       vtkm::Id& cellId,
                                       vtkm::Vec3f& parametric) const
  {
    return locator.FindCell(point, cellId, parametric);
  }
};

struct FindCellFunctor2
{
  template <typename Locator>
  VTKM_EXEC vtkm::ErrorCode operator()(Locator&& locator,
                                       const vtkm::Vec3f& point,
                                       vtkm::Id& cellId,
                                       vtkm::Vec3f& parametric,
                                       vtkm::Id3& prevCell) const
  {
    return locator.FindCell(point, cellId, parametric, prevCell);
  }
};

struct GetLeafDimensions
{
  template <typename Locator>
  VTKM_EXEC ReadPortal<DimVec3> operator()(Locator&& locator) const
  {
    return locator.GetLeafDimensions();
  }
};

struct GetLeafStartIndex
{
  template <typename Locator>
  VTKM_EXEC ReadPortal<vtkm::Id> operator()(Locator&& locator) const
  {
    return locator.GetLeafStartIndex();
  }
};

struct GetCellStartIndex
{
  template <typename Locator>
  VTKM_EXEC ReadPortal<vtkm::Id> operator()(Locator&& locator) const
  {
    return locator.GetCellStartIndex();
  }
};

struct GetCellCount
{
  template <typename Locator>
  VTKM_EXEC ReadPortal<vtkm::Id> operator()(Locator&& locator) const
  {
    return locator.GetCellCount();
  }
};

struct GetCellIds
{
  template <typename Locator>
  VTKM_EXEC ReadPortal<vtkm::Id> operator()(Locator&& locator) const
  {
    return locator.GetCellIds();
  }
};

struct GetTopLevel
{
  template <typename Locator>
  VTKM_EXEC vtkm::exec::Grid
  operator()(Locator&& locator) const
  {
    return locator.GetTopLevel();
  }
};
/*
template <typename CoordsPortalType>
struct GetCoords
{
  template <typename Locator>
  VTKM_EXEC CoordsPortalType
  operator()(Locator&& locator) const
  {
    return locator.GetCoords();
  }
};
*/

} // namespace detail

template <typename... LocatorTypes>
class VTKM_ALWAYS_EXPORT CellLocatorMultiplexer
{
using DimensionType = vtkm::Int16;
using DimVec3 = vtkm::Vec<DimensionType, 3>;
using FloatVec3 = vtkm::Vec3f;

template <typename T>
using ReadPortal = typename vtkm::cont::ArrayHandle<T>::ReadPortalType;

public:
  vtkm::exec::internal::Variant<LocatorTypes...> Locators;

public:
  CellLocatorMultiplexer() = default;

  template <typename Locator>
  VTKM_CONT CellLocatorMultiplexer(const Locator& locator)
    : Locators(locator)
  {
  }

#if 0
  VTKM_EXEC vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                                     vtkm::Id& cellId,
                                     vtkm::Vec3f& parametric,
                                     vtkm::Id3& /*prevCell*/) const
  {
    return this->FindCell(point, cellId, parametric);
  }
#endif

  VTKM_EXEC vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                                     vtkm::Id& cellId,
                                     vtkm::Vec3f& parametric) const
  {
    return this->Locators.CastAndCall(detail::FindCellFunctor{}, point, cellId, parametric);
  }

  VTKM_EXEC vtkm::ErrorCode FindCell(const vtkm::Vec3f& point,
                                     vtkm::Id& cellId,
                                     vtkm::Vec3f& parametric,
                                     vtkm::Id3& prevCell) const
  {
    return this->Locators.CastAndCall(detail::FindCellFunctor2{}, point, cellId, parametric, prevCell);
  }

  VTKM_EXEC ReadPortal<DimVec3> GetLeafDimensions() const { return this->Locators.CastAndCall(detail::GetLeafDimensions{}); }
  VTKM_EXEC ReadPortal<vtkm::Id> GetLeafStartIndex() const { return this->Locators.CastAndCall(detail::GetLeafStartIndex{}); }
  VTKM_EXEC ReadPortal<vtkm::Id> GetCellStartIndex() const { return this->Locators.CastAndCall(detail::GetCellStartIndex{}); }
  VTKM_EXEC ReadPortal<vtkm::Id> GetCellCount() const { return this->Locators.CastAndCall(detail::GetCellCount{}); }
  VTKM_EXEC ReadPortal<vtkm::Id> GetCellIds() const { return this->Locators.CastAndCall(detail::GetCellIds{}); }
  VTKM_EXEC vtkm::exec::Grid GetTopLevel() const {return this->Locators.CastAndCall(detail::GetTopLevel{});}

//  template <typename CoordsPortalType>
//  VTKM_EXEC CoordsPortalType GetCoords() const {return this->Locators.CastAndCall(detail::GetCoords{});

  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC CellLocatorMultiplexer* operator->() { return this; }
  VTKM_DEPRECATED(1.6, "Locators are no longer pointers. Use . operator.")
  VTKM_EXEC const CellLocatorMultiplexer* operator->() const { return this; }
};

}
} // namespace vtkm::exec

#endif //vtk_m_exec_CellLocatorMultiplexer_h

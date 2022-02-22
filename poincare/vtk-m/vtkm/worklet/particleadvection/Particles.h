//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_worklet_particleadvection_Particles_h
#define vtk_m_worklet_particleadvection_Particles_h

#include <vtkm/Geometry.h>
#include <vtkm/Particle.h>
#include <vtkm/Types.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ExecutionObjectBase.h>
#include <vtkm/worklet/particleadvection/IntegratorStatus.h>

namespace vtkm
{
namespace worklet
{
namespace particleadvection
{
template <typename ParticleType>
class ParticleExecutionObject
{
public:
  VTKM_EXEC_CONT
  ParticleExecutionObject()
    : Particles()
    , MaxSteps(0)
  {
  }

  ParticleExecutionObject(vtkm::cont::ArrayHandle<ParticleType> particleArray,
                          vtkm::Id maxSteps,
                          vtkm::cont::DeviceAdapterId device,
                          vtkm::cont::Token& token)
  {
    this->Particles = particleArray.PrepareForInPlace(device, token);
    this->MaxSteps = maxSteps;
  }

  VTKM_EXEC
  ParticleType GetParticle(const vtkm::Id& idx) { return this->Particles.Get(idx); }

  VTKM_EXEC
  void PreStepUpdate(const vtkm::Id& vtkmNotUsed(idx)) {}

  VTKM_EXEC
  void StepUpdate(const vtkm::Id& idx, vtkm::FloatDefault time, const vtkm::Vec3f& pt)
  {
    ParticleType p = this->GetParticle(idx);
    p.Pos = pt;
    p.Time = time;
    p.NumSteps++;
    this->Particles.Set(idx, p);
  }

  VTKM_EXEC
  void StatusUpdate(const vtkm::Id& idx,
                    const vtkm::worklet::particleadvection::IntegratorStatus& status,
                    vtkm::Id maxSteps)
  {
    ParticleType p = this->GetParticle(idx);

    if (p.NumSteps == maxSteps)
      p.Status.SetTerminate();

    if (status.CheckFail())
      p.Status.SetFail();
    if (status.CheckSpatialBounds())
      p.Status.SetSpatialBounds();
    if (status.CheckTemporalBounds())
      p.Status.SetTemporalBounds();
    if (status.CheckInGhostCell())
      p.Status.SetInGhostCell();
    this->Particles.Set(idx, p);
  }

  VTKM_EXEC
  bool CanContinue(const vtkm::Id& idx)
  {
    ParticleType p = this->GetParticle(idx);

    return (p.Status.CheckOk() && !p.Status.CheckTerminate() && !p.Status.CheckSpatialBounds() &&
            !p.Status.CheckTemporalBounds() && !p.Status.CheckInGhostCell());
  }

  VTKM_EXEC
  void UpdateTookSteps(const vtkm::Id& idx, bool val)
  {
    ParticleType p = this->GetParticle(idx);
    if (val)
      p.Status.SetTookAnySteps();
    else
      p.Status.ClearTookAnySteps();
    this->Particles.Set(idx, p);
  }

protected:
  using ParticlePortal = typename vtkm::cont::ArrayHandle<ParticleType>::WritePortalType;

  ParticlePortal Particles;
  vtkm::Id MaxSteps;
};

template <typename ParticleType>
class Particles : public vtkm::cont::ExecutionObjectBase
{
public:
  VTKM_CONT vtkm::worklet::particleadvection::ParticleExecutionObject<ParticleType>
  PrepareForExecution(vtkm::cont::DeviceAdapterId device, vtkm::cont::Token& token) const
  {
    return vtkm::worklet::particleadvection::ParticleExecutionObject<ParticleType>(
      this->ParticleArray, this->MaxSteps, device, token);
  }

  VTKM_CONT
  Particles(vtkm::cont::ArrayHandle<ParticleType>& pArray, vtkm::Id& maxSteps)
    : ParticleArray(pArray)
    , MaxSteps(maxSteps)
  {
  }

  Particles() {}

protected:
  vtkm::cont::ArrayHandle<ParticleType> ParticleArray;
  vtkm::Id MaxSteps;
};


template <typename ParticleType>
class StateRecordingParticleExecutionObject : public ParticleExecutionObject<ParticleType>
{
public:
  VTKM_EXEC_CONT
  StateRecordingParticleExecutionObject()
    : ParticleExecutionObject<ParticleType>()
    , History()
    , Length(0)
    , StepCount()
    , ValidPoint()
  {
  }

  StateRecordingParticleExecutionObject(vtkm::cont::ArrayHandle<ParticleType> pArray,
                                        vtkm::cont::ArrayHandle<vtkm::Vec3f> historyArray,
                                        vtkm::cont::ArrayHandle<vtkm::Id> validPointArray,
                                        vtkm::cont::ArrayHandle<vtkm::Id> stepCountArray,
                                        vtkm::Id maxSteps,
                                        vtkm::cont::DeviceAdapterId device,
                                        vtkm::cont::Token& token)
    : ParticleExecutionObject<ParticleType>(pArray, maxSteps, device, token)
    , Length(maxSteps + 1)
  {
    vtkm::Id numPos = pArray.GetNumberOfValues();
    this->History = historyArray.PrepareForOutput(numPos * this->Length, device, token);
    this->ValidPoint = validPointArray.PrepareForInPlace(device, token);
    this->StepCount = stepCountArray.PrepareForInPlace(device, token);
  }

  VTKM_EXEC
  void PreStepUpdate(const vtkm::Id& idx)
  {
    if (this->StepCount.Get(idx) == 0)
    {
      ParticleType p = this->ParticleExecutionObject<ParticleType>::GetParticle(idx);
      vtkm::Id loc = idx * this->Length;
      this->History.Set(loc, p.Pos);
      this->ValidPoint.Set(loc, 1);
      this->StepCount.Set(idx, 1);
    }
  }

  VTKM_EXEC
  void StepUpdate(const vtkm::Id& idx, vtkm::FloatDefault time, const vtkm::Vec3f& pt)
  {
    this->ParticleExecutionObject<ParticleType>::StepUpdate(idx, time, pt);

    //local step count.
    vtkm::Id stepCount = this->StepCount.Get(idx);

    vtkm::Id loc = idx * this->Length + stepCount;
    this->History.Set(loc, pt);
    this->ValidPoint.Set(loc, 1);
    this->StepCount.Set(idx, stepCount + 1);
  }

protected:
  using IdPortal = typename vtkm::cont::ArrayHandle<vtkm::Id>::WritePortalType;
  using HistoryPortal = typename vtkm::cont::ArrayHandle<vtkm::Vec3f>::WritePortalType;

  HistoryPortal History;
  vtkm::Id Length;
  IdPortal StepCount;
  IdPortal ValidPoint;
};

template <typename ParticleType>
class StateRecordingParticles : vtkm::cont::ExecutionObjectBase
{
public:
  //Helper functor for compacting history
  struct IsOne
  {
    template <typename T>
    VTKM_EXEC_CONT bool operator()(const T& x) const
    {
      return x == T(1);
    }
  };


  VTKM_CONT vtkm::worklet::particleadvection::StateRecordingParticleExecutionObject<ParticleType>
  PrepareForExecution(vtkm::cont::DeviceAdapterId device, vtkm::cont::Token& token) const
  {
    return vtkm::worklet::particleadvection::StateRecordingParticleExecutionObject<ParticleType>(
      this->ParticleArray,
      this->HistoryArray,
      this->ValidPointArray,
      this->StepCountArray,
      this->MaxSteps,
      device,
      token);
  }
  VTKM_CONT
  StateRecordingParticles(vtkm::cont::ArrayHandle<ParticleType>& pArray, const vtkm::Id& maxSteps)
    : MaxSteps(maxSteps)
    , ParticleArray(pArray)
  {
    vtkm::Id numParticles = static_cast<vtkm::Id>(pArray.GetNumberOfValues());

    //Create ValidPointArray initialized to zero.
    vtkm::cont::ArrayHandleConstant<vtkm::Id> tmp(0, (this->MaxSteps + 1) * numParticles);
    vtkm::cont::ArrayCopy(tmp, this->ValidPointArray);

    //Create StepCountArray initialized to zero.
    vtkm::cont::ArrayHandleConstant<vtkm::Id> tmp2(0, numParticles);
    vtkm::cont::ArrayCopy(tmp2, this->StepCountArray);
  }

  VTKM_CONT
  StateRecordingParticles(vtkm::cont::ArrayHandle<ParticleType>& pArray,
                          vtkm::cont::ArrayHandle<vtkm::Vec3f>& historyArray,
                          vtkm::cont::ArrayHandle<vtkm::Id>& validPointArray,
                          vtkm::Id& maxSteps)
  {
    this->ParticleArray = pArray;
    this->HistoryArray = historyArray;
    this->ValidPointArray = validPointArray;
    this->MaxSteps = maxSteps;
  }

  VTKM_CONT
  void GetCompactedHistory(vtkm::cont::ArrayHandle<vtkm::Vec3f>& positions)
  {
    vtkm::cont::Algorithm::CopyIf(this->HistoryArray, this->ValidPointArray, positions, IsOne());
  }

protected:
  vtkm::cont::ArrayHandle<vtkm::Vec3f> HistoryArray;
  vtkm::Id MaxSteps;
  vtkm::cont::ArrayHandle<ParticleType> ParticleArray;
  vtkm::cont::ArrayHandle<vtkm::Id> StepCountArray;
  vtkm::cont::ArrayHandle<vtkm::Id> ValidPointArray;
};


template <typename ParticleType>
class PoincareParticleExecutionObject : public ParticleExecutionObject<ParticleType>
{
public:
  VTKM_EXEC_CONT
  PoincareParticleExecutionObject()
    : ParticleExecutionObject<ParticleType>()
    , History()
    , Length(0)
    , MaxPunctures(0)
    , Plane()
    , PunctureCount()
    , StepCount()
    , ValidPoint()
  {
  }

  PoincareParticleExecutionObject(vtkm::cont::ArrayHandle<ParticleType> pArray,
                                  vtkm::cont::ArrayHandle<vtkm::Vec3f> historyArray,
                                  vtkm::cont::ArrayHandle<vtkm::Id> validPointArray,
                                  vtkm::cont::ArrayHandle<vtkm::Id> stepCountArray,
                                  vtkm::cont::ArrayHandle<vtkm::Id> punctureCountArray,
                                  vtkm::Id maxSteps,
                                  vtkm::Id maxPunctures,
                                  const vtkm::Plane<>& plane,
                                  vtkm::cont::DeviceAdapterId device,
                                  vtkm::cont::Token& token)
    : ParticleExecutionObject<ParticleType>(pArray, maxSteps, device, token)
    , Length(maxPunctures + 1)
    , MaxPunctures(maxPunctures)
    , Plane(plane)
  {
    vtkm::Id numPos = pArray.GetNumberOfValues();
    this->History = historyArray.PrepareForOutput(numPos * this->Length, device, token);
    this->ValidPoint = validPointArray.PrepareForInPlace(device, token);
    this->StepCount = stepCountArray.PrepareForInPlace(device, token);
    this->PunctureCount = punctureCountArray.PrepareForInPlace(device, token);
  }

  VTKM_EXEC
  void StatusUpdate(const vtkm::Id& idx,
                    const vtkm::worklet::particleadvection::IntegratorStatus& status,
                    vtkm::Id maxSteps)
  {
    ParticleType p = this->GetParticle(idx);
    vtkm::Id puncCount = this->PunctureCount.Get(idx);

    if (p.NumSteps == maxSteps)
      p.Status.SetTerminate();
    if (puncCount >= this->MaxPunctures)
      p.Status.SetTerminate();

    if (status.CheckFail())
      p.Status.SetFail();
    if (status.CheckSpatialBounds())
      p.Status.SetSpatialBounds();
    if (status.CheckTemporalBounds())
      p.Status.SetTemporalBounds();
    if (status.CheckInGhostCell())
      p.Status.SetInGhostCell();
    this->Particles.Set(idx, p);
  }

  VTKM_EXEC
  void PreStepUpdate(const vtkm::Id& idx)
  {
    if (this->StepCount.Get(idx) == 0)
    {
      ParticleType p = this->ParticleExecutionObject<ParticleType>::GetParticle(idx);
      vtkm::Id loc = idx * this->Length;
      this->History.Set(loc, p.Pos);
      this->ValidPoint.Set(loc, 1);
      this->StepCount.Set(idx, 1);
    }
  }

  VTKM_EXEC
  void StepUpdate(const vtkm::Id& idx, vtkm::FloatDefault time, const vtkm::Vec3f& pt)
  {
    //Get the previous point and calculate the distance to the plane.
    ParticleType prevPt = this->GetParticle(idx);
    auto prevDist = this->Plane.DistanceTo(prevPt.Pos);
    auto currDist = this->Plane.DistanceTo(pt);
    this->ParticleExecutionObject<ParticleType>::StepUpdate(idx, time, pt);

    //std::cout<<"  StepUpdate: "<<prevPt.Pos<<" --> "<<pt<<"  :: "<<prevDist<<" "<<currDist<<std::endl;

    //local step count.
    vtkm::Id stepCount = this->StepCount.Get(idx);

    //If we hit the plane, update the history.
    auto prevBit = vtkm::SignBit(prevDist), currBit = vtkm::SignBit(currDist);
    if (!prevBit && currBit)
    {
      vtkm::Id puncCount = this->PunctureCount.Get(idx);
      vtkm::Id loc = idx * this->Length + puncCount;
      this->History.Set(loc, pt);
      this->ValidPoint.Set(loc, 1);
      this->PunctureCount.Set(idx, puncCount + 1);
      //std::cout<<"  "<<idx<<": Poinc::StepUpdate: "<<puncCount<<" : "<<stepCount<<" "<<prevPt.Pos<<" ==> "<<pt<<" d= "<<prevDist<<" "<<currDist<<std::endl;
    }

    this->StepCount.Set(idx, stepCount + 1);

  }

protected:
  using IdPortal = typename vtkm::cont::ArrayHandle<vtkm::Id>::WritePortalType;
  using HistoryPortal = typename vtkm::cont::ArrayHandle<vtkm::Vec3f>::WritePortalType;

  HistoryPortal History;
  vtkm::Id Length;
  vtkm::Id MaxPunctures;
  vtkm::Plane<> Plane;
  IdPortal PunctureCount;
  IdPortal StepCount;
  IdPortal ValidPoint;
};


template <typename ParticleType>
class PoincareParticles : vtkm::cont::ExecutionObjectBase
{
public:
  //Helper functor for compacting history
  struct IsOne
  {
    template <typename T>
    VTKM_EXEC_CONT bool operator()(const T& x) const
    {
      return x == T(1);
    }
  };


  VTKM_CONT vtkm::worklet::particleadvection::PoincareParticleExecutionObject<ParticleType>
  PrepareForExecution(vtkm::cont::DeviceAdapterId device, vtkm::cont::Token& token) const
  {
    return vtkm::worklet::particleadvection::PoincareParticleExecutionObject<ParticleType>(
      this->ParticleArray,
      this->HistoryArray,
      this->ValidPointArray,
      this->StepCountArray,
      this->PunctureCountArray,
      this->MaxSteps,
      this->MaxPunctures,
      this->Plane,
      device,
      token);
  }
  VTKM_CONT
  PoincareParticles(vtkm::cont::ArrayHandle<ParticleType>& pArray, const vtkm::Plane<>& plane, const vtkm::Id& maxSteps, const vtkm::Id& maxPunctures)
    : MaxSteps(maxSteps)
    , MaxPunctures(maxPunctures)
    , ParticleArray(pArray)
    , Plane(plane)
  {
    vtkm::Id numParticles = static_cast<vtkm::Id>(pArray.GetNumberOfValues());

    //Create ValidPointArray initialized to zero.
    vtkm::cont::ArrayHandleConstant<vtkm::Id> tmp(0, (this->MaxPunctures + 1) * numParticles);
    vtkm::cont::ArrayCopy(tmp, this->ValidPointArray);

    //Create StepCountArray initialized to zero.
    vtkm::cont::ArrayHandleConstant<vtkm::Id> tmp2(0, numParticles);
    vtkm::cont::ArrayCopy(tmp2, this->StepCountArray);

    //Create PunctureCountArray initialized to zero.
    vtkm::cont::ArrayCopy(tmp2, this->PunctureCountArray);
  }

  VTKM_CONT
  PoincareParticles(vtkm::cont::ArrayHandle<ParticleType>& pArray,
                    vtkm::cont::ArrayHandle<vtkm::Vec3f>& historyArray,
                    vtkm::cont::ArrayHandle<vtkm::Id>& validPointArray,
                    const vtkm::Plane<>& plane,
                    vtkm::Id& maxSteps,
                    const vtkm::Id& maxPunctures)
    : Plane(plane)
  {
    this->ParticleArray = pArray;
    this->HistoryArray = historyArray;
    this->ValidPointArray = validPointArray;
    this->MaxSteps = maxSteps;
    this->MaxPunctures = maxPunctures;
  }

  VTKM_CONT
  void GetCompactedHistory(vtkm::cont::ArrayHandle<vtkm::Vec3f>& positions)
  {
    vtkm::cont::Algorithm::CopyIf(this->HistoryArray, this->ValidPointArray, positions, IsOne());
  }

  VTKM_CONT
  vtkm::cont::ArrayHandle<vtkm::Id> GetPunctureCountArray() const
  {
    return this->PunctureCountArray;
  }

protected:
  vtkm::cont::ArrayHandle<vtkm::Vec3f> HistoryArray;
  vtkm::Id MaxSteps;
  vtkm::Id MaxPunctures;
  vtkm::cont::ArrayHandle<ParticleType> ParticleArray;
  vtkm::cont::ArrayHandle<vtkm::Id> StepCountArray;
  vtkm::cont::ArrayHandle<vtkm::Id> PunctureCountArray;
  vtkm::cont::ArrayHandle<vtkm::Id> ValidPointArray;
  vtkm::Plane<> Plane;
};


} //namespace particleadvection
} //namespace worklet
} //namespace vtkm

#endif // vtk_m_worklet_particleadvection_Particles_h
//============================================================================

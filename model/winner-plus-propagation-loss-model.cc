/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011, 2012 Centre Tecnologic de Telecomunicacions de Catalunya (CTTC)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Marco Miozzo  <marco.miozzo@cttc.es>,
 *         Nicola Baldo <nbaldo@cttc.es>beep
 * 
 */
#include "ns3/log.h"
#include "ns3/double.h"
#include "ns3/boolean.h"
#include "ns3/enum.h"
#include "ns3/mobility-model.h"
#include <cmath>

#include "winner-plus-propagation-loss-model.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("WinnerPlusPropagationLossModel");

NS_OBJECT_ENSURE_REGISTERED (WinnerPlusPropagationLossModel);


TypeId
WinnerPlusPropagationLossModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::WinnerPlusPropagationLossModel")
    .SetParent<PropagationLossModel> ()
    .SetGroupName ("Propagation")
    .AddConstructor<WinnerPlusPropagationLossModel> ()
    .AddAttribute ("Frequency",
                   "The propagation frequency in Hz",
                   DoubleValue (806e6),
                   MakeDoubleAccessor (&WinnerPlusPropagationLossModel::m_frequency),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("LineOfSight",
                   "Whether or not line of sight is available",
                   BooleanValue (false),
                   MakeBooleanAccessor (&WinnerPlusPropagationLossModel::m_los),
                   MakeBooleanChecker ())
    .AddAttribute ("Environment",
                   "Environment Scenario",
                   EnumValue (UMaEnvironment),
                   MakeEnumAccessor (&WinnerPlusPropagationLossModel::m_environment),
                   MakeEnumChecker (UMiEnvironment, "UrbanMicro",
                                    O2IaEnvironment, "OutdoorToIndoorA",
                                    SMaEnvironment, "SuburbanMacro",
                                    UMaEnvironment, "UrbanMacro",
                                    O2IbEnvironment, "OutdoorToIndoorB"))
    .AddAttribute ("Layout",
                   "Grid Layout",
                   EnumValue (Hexagonal),
                   MakeEnumAccessor (&WinnerPlusPropagationLossModel::m_layout),
                   MakeEnumChecker (Hexagonal, "HexagonalLayout")); // Manhattan TBD
                                    
  return tid;
}

WinnerPlusPropagationLossModel::WinnerPlusPropagationLossModel ()
  : PropagationLossModel ()
{
}

WinnerPlusPropagationLossModel::~WinnerPlusPropagationLossModel ()
{
}

double
WinnerPlusPropagationLossModel::GetLoss (Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  double loss = 0.0;
  double fghz = m_frequency / 1e9;
  NS_ASSERT_MSG (fghz >= 0.45 && fghz <= 6, "Frequency must be between 0.45 GHz and 6.0 GHz");
  double dist = a->GetDistanceFrom (b) / 1000.0; 
  // double hbs = (a->GetPosition ().z > b->GetPosition ().z ? a->GetPosition ().z : b->GetPosition ().z);
  double hbs = 10; // Defined in WinnerPlus D5.3 v 1.0 Table 4-1
  // double hms = (a->GetPosition ().z < b->GetPosition ().z ? a->GetPosition ().z : b->GetPosition ().z);
  double hms = 1.5; // Defined in WinnerPlus D5.3 v 1.0 Table 4-1

  double dbp = 4*hbs*hMS*fghz/3e8;
  double h_tick_bs = hbs - 1.0;
  double h_tick_ms = hms - 1.0;
  double d_tick_bp = 4*h_tick_bs*h_tick_MS*fghz/3e8;
  if (m_environment == UMiEnvironment){
    if (m_los == true){
      loss = 22.7 * std::log10(dist) + 27.0 + 20.0 * std::log10(fghz);
    }
    else if (m_los == false){
      // Currently only hexagonal Layout is supported
      if (fghz < 1.5){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 16.33 + 26.16 * std::log10(fghz);
      }
      else if (fghz < 2.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 14.78 + 34.97 * std::log10(fghz);
      }
      else if (fghz < 6.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 18.38 + 23.00 * std::log10(fghz);
      }
    }
  }
  else if (m_environment == O2IaEnvironment){
    double din = 0.003; // Let's assume that the indoor range is 3m
    double dout = dist - din;
    if (fghz < 1.5){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 16.33 + 26.16 * std::log10(fghz);
    }
    else if (fghz < 2.0){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 14.78 + 34.97 * std::log10(fghz);
    }
    else if (fghz < 6.0){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 18.38 + 23.00 * std::log10(fghz);
    }
    loss = PLb*(dout + din) + 21.04 + 14*(1-1.8*std::log10(fghz)) + 0.5*din - 0.8*hms;
  }
  else if (m_environment == SMaEnvironment){
    if (m_los == true){
      if (dist < dbp){
        loss = 23.8 * std::log10(dist) + 27.2 + 20.0 * std::log10(fghz);
      }
      else{
        loss = 40.0 * std::log10(dist) + 9.0 - 16.2*std::log10(hbs) -16.2*std::log10(hms) + 3.8*std::log10(fghz);
      }
    }
    else if (m_los == false){
      if (fghz < 1.5){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 13.33 + 26.16 * std::log10(fghz);
      }
      else if (fghz < 2.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 11.78 + 34.97 * std::log10(fghz);
      }
      else if (fghz < 6.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 15.38 + 23.00 * std::log10(fghz);
      }
    }
  }
  else if (m_environment == UMaEnvironment){
    if (m_los == true){
      if (dist < d_tick_bp){
        loss = 26.0 * std::log10(dist) + 25.0 + 20.0 * std::log10(fghz);
      }
      else{
        loss = 40.0 * std::log10(dist) + 9.27 - 14.0*std::log10(hbs) -14.0*std::log10(hms) + 6.0*std::log10(fghz);
      }
    }
    else if (m_los == false){
      if (fghz < 1.5){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 16.33 + 26.16 * std::log10(fghz);
      }
      else if (fghz < 2.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 14.78 + 34.97 * std::log10(fghz);
      }
      else if (fghz < 6.0){
        loss = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 18.38 + 23.00 * std::log10(fghz);
      }
    }
  }
  else if (m_environment == O2IbEnvironment){
    double din = 0.003 // Let's assume that the indoor range is 3m
    double dout = dist - din
    if (fghz < 1.5){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 16.33 + 26.16 * std::log10(fghz);
    }
    else if (fghz < 2.0){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 14.78 + 34.97 * std::log10(fghz);
    }
    else if (fghz < 6.0){
      double PLb = (44.9 - 6.55 * std::log10(hbs)) * std::log10(dist) + 5.83 * std::log10(hbs) + 18.38 + 23.00 * std::log10(fghz);
    }
    loss = PLb*(dout + din)) + 21.04 + 14*(1-1.8*std::log10(fghz)) + 0.5*din - 0.8*hms;
  }
  return loss;
}

double 
WinnerPlusPropagationLossModel::DoCalcRxPower (double txPowerDbm,
						Ptr<MobilityModel> a,
						Ptr<MobilityModel> b) const
{
  return (txPowerDbm - GetLoss (a, b));
}

int64_t
WinnerPlusPropagationLossModel::DoAssignStreams (int64_t stream)
{
  return 0;
}


} // namespace ns3

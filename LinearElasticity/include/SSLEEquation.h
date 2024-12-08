#pragma once

#include <talyfem/fem/cequation.h>
#include <talyfem/stabilizer/tezduyar_upwind.h>
#include "LENodeData.h"
#include "LEInputData.h"
#include <talyfem/talyfem.h>
#include <link.h>
#include <cmath>
#include "SBMCalc.h"

class SSLEEquation : public TALYFEMLIB::CEquation<LENodeData>
{

  const bool IFSBM = true;
  bool ShortestDist = true;
  const IBM_METHOD method;
  const IMGA *imga_;
  my_kd_tree_t *kd_tree_;


public:
  explicit SSLEEquation(LEInputData *idata, const IBM_METHOD _method, const IMGA *imga)
      : TALYFEMLIB::CEquation<LENodeData>(false, TALYFEMLIB::kAssembleGaussPoints), method(_method), imga_(imga)
  {
    idata_ = idata;
  }

  void Solve(double dt, double t) override
  {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override
  {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae)
  {

    double thickness = 1;

    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = fe.nsd();

    // # of basis functions
    const int n_basis_functions = fe.nbf();

    const double detJxW = fe.detJxW();

#if (DIM == 2)

    //    const double detJxW = fe.detJxW();
    for (int a = 0; a < n_basis_functions; a++)
    {
      for (int b = 0; b < n_basis_functions; b++)
      {
        for (int i = 0; i < DIM; i++)
        {
          for (int j = 0; j < DIM; j++)
          {
            Ae(DIM * a + i, DIM * b + j) += fe.dN(b, j) * idata_->Cmatrix[i][j] * fe.dN(a, i) * detJxW + fe.dN(b, DIM - j - 1) * idata_->Cmatrix[DIM][DIM] * fe.dN(a, DIM - i - 1) * detJxW;
          }
        }
      }
    }

#endif

#if (DIM == 3)

    //            std::vector<std::vector<double>> Ae_check(n_dimensions * n_basis_functions, std::vector<double>(n_dimensions * n_basis_functions));
    //            std::vector<std::vector<double>> Be(n_dimensions * n_basis_functions, std::vector<double>(6));
    //            CalcBe(fe,Be);
    //
    //            std::vector<std::vector<double>> BeCmatrix(n_dimensions * n_basis_functions, std::vector<double>(6));
    //            CalcBeCmatrix(fe, Be, idata_->Cmatrix, BeCmatrix);
    //
    ////             for final term -> previous implementation!
    //            for (int a = 0; a < n_dimensions * n_basis_functions; a++)
    //            {
    //                for (int b = 0; b < n_dimensions * n_basis_functions; b++)
    //                {
    ////                    Ae_check[a][b] = 0;
    //                    for (int k = 0; k < 6; k++)
    //                    {
    //                        Ae(a, b) += BeCmatrix[a][k] * Be[b][k] * detJxW;
    ////                        Ae_check[a][b] += BeCmatrix[a][k] * Be[b][k] * detJxW;
    //                    }
    //                }
    //            }

    /*
     * below is from CalcAe3D.m.
     * And, we try to know the pattern by ourselves
     */
    const double Emv = idata_->Cmatrix[0][0];
    const double Ev = idata_->Cmatrix[1][0];
    const double half = idata_->Cmatrix[3][3];
    //
    //    for (int i = 0; i < fe.nbf(); ++i) {
    //        for (int j = 0; j < fe.nbf(); ++j) {
    //            for (int k1 = 0; k1 < DIM; k1++) {
    //                for (int k2 = 0; k2 < DIM; k2++) {
    //                    // First type of element
    //                    if (k1 == k2) {
    //                        Ae(DIM * i + k1, DIM * j + k2) += Emv * fe.dN(i, k1) * fe.dN(j, k1) * detJxW;
    ////                        Ae_check[DIM * i + k1][ DIM * j + k2]-= Emv * fe.dN(i, k1) * fe.dN(j, k1) * detJxW;
    //                        if (k1 == 0) {
    //                            Ae(DIM * i + k1, DIM * j + k2) +=
    //                                    half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
    ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
    //                        } else if (k1 == 1) {
    //                            Ae(DIM * i + k1, DIM * j + k2) +=
    //                                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
    ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
    //                        } else if (k1 == 2) {
    //                            Ae(DIM * i + k1, DIM * j + k2) +=
    //                                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
    ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
    //                        }
    //                    }
    //                        // Second type of element -> switch k1 and k2
    //                    else {
    //                        Ae(DIM * i + k1, DIM * j + k2) +=
    //                                (Ev * fe.dN(i, k1) * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
    ////                        Ae_check[DIM * i + k1][ DIM * j + k2]-= (Ev * fe.dN(i, k1) * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
    //                    }
    //                }
    //            }
    //        }
    //    }

    const int nbf = fe.nbf();
    for (int i = 0; i < nbf; ++i)
    {
      int DIM_i = DIM * i;

      for (int j = 0; j < nbf; ++j)
      {
        int DIM_j = DIM * j;

        for (int k1 = 0; k1 < DIM; k1++)
        {
          double dN_ik1 = fe.dN(i, k1); // Cache this value for k1 loop

          for (int k2 = 0; k2 < DIM; k2++)
          {
            if (k1 == k2)
            {
              Ae(DIM_i + k1, DIM_j + k2) += Emv * dN_ik1 * fe.dN(j, k1) * detJxW;

              switch (k1)
              {
              case 0:
                Ae(DIM_i + k1, DIM_j + k2) +=
                    half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
                break;
              case 1:
                Ae(DIM_i + k1, DIM_j + k2) +=
                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
                break;
              case 2:
                Ae(DIM_i + k1, DIM_j + k2) +=
                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
                break;
              }
            }
            else
            {
              Ae(DIM_i + k1, DIM_j + k2) +=
                  (Ev * dN_ik1 * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
            }
          }
        }
      }
    }

#endif
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be)
  {
    return;
    using namespace TALYFEMLIB;
    const ZEROPTV &p = fe.position();

    ZEROPTV BodyForce;
    bool ForceHaveSet = false;
    CalcForce(p, BodyForce, ForceHaveSet);
    if (!ForceHaveSet)
    {
      BodyForce = idata_->BodyForce;
    }

#if (DIM == 3)
    double body_z = idata_->BodyForce[2];
#endif
    double BR_V = idata_->radialbodyforce.br_v;
    int BR_POW = idata_->radialbodyforce.br_pow;

//      BodyForce.print();

    /*
     * please write something like this so that it can support both 2D and 3D
     */
    for (int a = 0; a < fe.nbf(); a++)
    {
      for (int dim = 0; dim < DIM; dim++)
      {
        be(DIM * a + dim) += fe.N(a) * BodyForce(dim) * fe.detJxW();
      }
    }

#if (DIM == 2)
    double x_min = idata_->mesh_def.physDomain.min[0];
    double y_min = idata_->mesh_def.physDomain.min[1];
    double x_max = idata_->mesh_def.physDomain.max[0];
    double y_max = idata_->mesh_def.physDomain.max[1];
    double x_mid = (x_min + x_max) / 2;
    double y_mid = (y_min + y_max) / 2;

    double x = p.x();
    double y = p.y();

    for (int a = 0; a < fe.nbf(); a++)
    {
      double r_value = 0;
      double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
      r_value = BR_V * pow(radius, BR_POW);
      // std::cout<<"r_value:"<<r_value<<"\n";
      double sin_x = (x - x_mid) / radius;
      double sin_y = (y - y_mid) / radius;
      be(2 * a) += fe.N(a) * r_value * sin_x * fe.detJxW();
      be(2 * a + 1) += fe.N(a) * r_value * sin_y * fe.detJxW();
    }
#endif
  }

  ///// ==================== ibm start====================================
#pragma mark Ae-be Surface IBM

  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {
    assert(method == IBM_METHOD::NITSCHE);
  } // end:ibm_Ae

  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {
    assert(method == IBM_METHOD::NITSCHE);
  } // end:ibm_be

  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                              const NodeAndValues<DENDRITE_REAL> &gpinfo,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h,
                              const std::vector<double> &surface_values)
  {
    ibm_Integrands4side_Ae(fe, Ae, gpinfo, position, h);
  }

  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &Ae,
                              const NodeAndValues<DENDRITE_REAL> &gpinfo,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h,
                              const std::vector<double> &surface_values)
  {
    ibm_Integrands4side_be(fe, Ae, gpinfo, position, h);
  }
  ///// ==================== ibm end====================================
#pragma mark Ae-be Surface Carved
  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZeroMatrix<double> &Ae) {

        assert(method == IBM_METHOD::SBM);


        double h = ElementSize(fe);
        const int n_basis_functions = fe.nbf();


        const double detSideJxW = fe.detJxW();

        double eps = 1e-15;
        bool x_minus_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.min[0]) < eps;
        bool y_minus_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.min[1]) < eps;
        bool x_plus_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.max[0]) < eps;
        bool y_plus_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.max[1]) < eps;

        if(idata_->NormalTraction.boundary_fitted)
        {
//            we apply nothing to the boundary fitted side in Ae matrix
//            we just return saying TSM
            return;
        }

        if(x_minus_wall || y_minus_wall || y_plus_wall)
        {
//        this is the case where the left boundary is fixed
            return;
        }
        double d[DIM];
        if(x_plus_wall)
        {
            d[0] = -h/idata_->elemOrder;

            d[1] = 0;

// THIS IS A SPECIAL CASE WHERE WE ARE APPLYING TRACTION TO THE TOP WALL
// if WE APPLY BC THROUGH RIGHT WALL WE
// WILL GET COMPETITION BETWEEN THE TWO BC
        }

        if(y_plus_wall)
        {
            d[0] = 0;
            d[1] = -h/idata_->elemOrder;
        }



        //////////////////////////////////////weak//////////////////////////////////////////

        double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];
        std::vector<std::vector<double>> Be(DIM * n_basis_functions);
//  gradStrainDotDistanceMatrix
        std::vector<std::vector<double>> gradStrainDotDistance(DIM * n_basis_functions);
        /// for middle term => B_T*C_T
        std::vector<std::vector<double>> BeCmatrix(DIM * n_basis_functions);
        double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
        /// for mid2 term => B_T*C_T*n
        std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);

        CalcCmatrix(Cmatrix);
//            Strain Dot Distance Matrix
        CalcGradStrainDotDistance(fe, gradStrainDotDistance, d);

//            multiply Cmatrix with gradStrainDotDistance
        std::vector<std::vector<double>> CMatrixGradStrainDotDistance(DIM * n_basis_functions);

//            just use THIS FUNCTION TO MULTIPLY THE TWO MAT IT'S ALREADY IMPLEMENTED
        CalcBeCmatrix(fe, gradStrainDotDistance, Cmatrix, CMatrixGradStrainDotDistance);
//            multiply CMatrixGradStrainDotDistance with normal
        CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);
//            multiply CMatrixGradStrainDotDistance with normal
        std::vector<std::vector<double>> Hessian_Dot_Normal(DIM * n_basis_functions);
//            Hessian_Dot_Normal is DIM * (n_basis_functions*DIM)
        CalcStressDotNormal(fe, CMatrixGradStrainDotDistance, SurrogateNormalMatrix, Hessian_Dot_Normal);


//            get the basis matrix , we are just doing this to make our life easier
        std::vector<std::vector<double>> basisMatrix(DIM * n_basis_functions);

        CalcBasisMatrix(fe, basisMatrix);
#ifndef NDEBUG
        for (int i = 0; i < DIM * n_basis_functions; i++) {
            for (int j = 0; j < DIM; j++) {
                if (std::isnan(basisMatrix[i][j])) {
                    std::cout << "basisMatrix is nan\n";
                }
            }
        }
        for (int i = 0; i < DIM * n_basis_functions; i++) {
          for (int j = 0; j < DIM; j++) {
              if (std::isnan(Hessian_Dot_Normal[i][j])) {
                  std::cout << "Hessian_Dot_Normal is nan\n";
              }
          }
      }
#endif

/// THIS IS THE MAIN INTEGRAL WE ARE LOOKING FOR


        for (int a = 0; a < DIM * n_basis_functions; a++) {
            for (int b = 0; b < DIM * n_basis_functions; b++) {
                double N = 0;
                for (int k = 0; k < DIM; k++) {
                    N += Hessian_Dot_Normal[a][k] * basisMatrix[b][k] * detSideJxW;


                }
                Ae(a, b) += N;
            }
        }

        return;



    }
  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZEROARRAY<double> &be)
  {
      double eps = 1e-10;
      double h = ElementSize(fe);
      bool x_minus_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.min[0]) < eps;
      bool y_minus_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.min[1]) < eps;
      bool x_plus_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.max[0]) < eps;
      bool y_plus_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.max[1]) < eps;
    if(x_minus_wall || y_minus_wall || y_plus_wall)
    {
//        this is the case where the left boundary is fixed
        return;
    }
    double d[DIM];

    if(x_plus_wall)
    {
        d[0] = -h;
        d[1] = 0;
//        THIS IS A SPECIAL CASE WHERE WE ARE APPLYING TRACTION TO THE TOP WALL
//        if WE APPLY BC THROUGH RIGHT WALL WE
//        WILL GET COMPETITION BETWEEN THE TWO BC
    }
    if(y_plus_wall)
    {
        d[0] = 0;
        d[1] = -h;
    }


  ZEROPTV traction = CalcNormalTraction(fe, d);


    for (int a = 0; a < fe.nbf(); a++)
    {
        for (int dim = 0; dim < DIM; dim++)
        {
            be(DIM * a + dim) += fe.N(a) * traction(dim) * fe.detJxW();
        }
    }

  }

#pragma mark Ae-be to support more dendrite-kt version

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h)
  {
    Integrands_Ae(fe, Ae);
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const double *h)
  {
    Integrands_be(fe, be);
  }

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae,
                          const double *h)
  {

    Integrands4side_Ae(fe, side_idx, id, Ae);
  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be,
                          const double *h)
  {
    Integrands4side_be(fe, side_idx, id, be);
  }

  void setSbmCalc(my_kd_tree_t *kd_tree)
    {
        kd_tree_ = kd_tree;
    }

private:
  LEInputData *idata_;

  /**
   * This function is going to check the SBMGeo and find out the corresponding MMS force calculated from MMS solution
   * @param p IN
   * @param BodyForce OUT
   * @param ForceHaveSet OUT
   */
  void CalcForce(const TALYFEMLIB::ZEROPTV &p, ZEROPTV &BodyForce, bool &ForceHaveSet) const
  {
    switch (idata_->SbmGeo)
    {

    case LEInputData::SBMGeo::CIRCLE:
    case LEInputData::SBMGeo::STAR:
    {
      ///
      double pi = M_PI;
      double E = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;

      ///
      double x = p.x();
      double y = p.y();

      ///
      BodyForce.x() = (E * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (5 * (2 * poisson + 2)) - (E * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1)) - (E * poisson * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1));
      BodyForce.y() = (E * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (5 * (2 * poisson + 2)) - (E * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1)) - (E * poisson * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1));
      ForceHaveSet = true;
      break;
    }

    case LEInputData::SBMGeo::ROTATE:
    {
      ///
      double pi = M_PI;
      double E = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;

      ///
      double x_rotate = p.x() - 0.5 - 0.21;
      double y_rotate = p.y() - 0.5 - 0.21;

      /// rotate back to [0,0] to [1,1]
      double x = cos(M_PI / 4) * x_rotate - sin(M_PI / 4) * y_rotate + 0.5;
      double y = sin(M_PI / 4) * x_rotate + cos(M_PI / 4) * y_rotate + 0.5;

      // std::cout<< "x,y = " << x << "," << y << "\n";

      ///
      double BodyForceXp = (E * pi * pi * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / 210 + (pi * pi * cos(pi * x) * sin(pi * y)) / 10)) / (2 * poisson + 2) + (E * poisson * pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / (210 * (pow(poisson, 2) - 1));
      double BodyForceYp = (E * poisson * pi * pi * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1)) - (E * pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / (90 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos(pi * y) * sin(pi * x)) / 10 - (pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / 490)) / (2 * poisson + 2);

      ///
      BodyForce.x() = (BodyForceXp * cos(M_PI / 4) + BodyForceYp * cos(M_PI / 4));
      BodyForce.y() = (BodyForceYp * sin(M_PI / 4) - BodyForceXp * sin(M_PI / 4));
      ForceHaveSet = true;
      break;
    }

#if (DIM == 3)
    case LEInputData::SBMGeo::SPHERE:
    {
      ///
      double pi = M_PI;
      double E = idata_->planeStress.young;
      double v = idata_->planeStress.poisson;

      ///
      double x = p.x();
      double y = p.y();
      double z = p.z();

      ///
      BodyForce.x() = (E * ((pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / 20 + (pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / 10) * (v - 0.5)) / ((2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / (20 * (2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z) * (v - 1)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z) * (v - 0.5)) / (5 * (2 * v - 1) * (v + 1));
      BodyForce.y() = (E * ((pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / 10 + (pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / 20) * (v - 0.5)) / ((2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / (10 * (2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / (20 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z) * (v - 1)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z) * (v - 0.5)) / (5 * (2 * v - 1) * (v + 1));
      BodyForce.z() = (E * v * pow(pi, 2) * cos(pi * x) * cos(pi * y) * cos(pi * z)) / (5 * (2 * v - 1) * (v + 1)) - (2 * E * (v - 0.5) * ((pow(pi, 2) * cos(pi * x) * cos(pi * y) * cos(pi * z)) / 10 - (pow(pi, 2) * cos(pi * z) * sin(pi * x) * sin(pi * y)) / 20)) / ((2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * z) * sin(pi * x) * sin(pi * y) * (v - 1)) / (20 * (2 * v - 1) * (v + 1));
      ForceHaveSet = true;
      break;
    }
#endif

    case LEInputData::SBMGeo::NONE:
    {
      ForceHaveSet = false;
      break;
    }

    default:
    {
      break;
    }
    }
  }

  DENDRITE_REAL normalDistance(const TALYFEMLIB::FEMElm &fe,
                               const ZEROPTV &normal,
                               const ZEROPTV &h)
  {
    double max_h = std::max(std::max(h.x(), h.y()), h.z());

    ZEROPTV root_node;
    fe.grid()->GetCoord(root_node, 0);
    int numNode = pow(idata_->elemOrder + 1, DIM);
    std::vector<ZEROPTV> elm_node(numNode);
#ifdef ENABLE_3D
    for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
    {
      for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
      {
        for (int p3 = 0; p3 < idata_->elemOrder + 1; p3++)
        {
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3] = root_node;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].x() += p3 * h.x() / idata_->elemOrder;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].y() += p2 * h.y() / idata_->elemOrder;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].z() += p1 * h.y() / idata_->elemOrder;
        }
      }
    }
#else
    for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
    {
      for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
      {
        elm_node[p1 * (idata_->elemOrder + 1) + p2] = root_node;
        elm_node[p1 * (idata_->elemOrder + 1) + p2].x() += p2 * h.x() / idata_->elemOrder;
        elm_node[p1 * (idata_->elemOrder + 1) + p2].y() += p1 * h.y() / idata_->elemOrder;
      }
    }
#endif
    std::vector<double> hb(numNode);
    for (unsigned int i = 0; i < numNode; i++)
    {
      ZEROPTV vec1 = fe.position() - elm_node[i];
      hb[i] = vec1.innerProduct(normal) > 0 ? vec1.innerProduct(normal) : 0;
    }
    auto maxhb = std::max_element(hb.begin(), hb.end());
    if (*maxhb < max_h * 1e-2)
    {
      return max_h * 1e-2;
    }
    return *maxhb;
  }

#pragma mark normal matrix

  /*
   * this normal matrix calculation is based on t_i = stress_ij * n_j
   * the total formulation please check cheng-hau's LE presentation
   */
  void CalcSurrogateNormalMatrix(const TALYFEMLIB::FEMElm &fe, double (&SurrogateNormalMatrix)[DIM][3 * (DIM - 1)])
  {
#if (DIM == 2)
    SurrogateNormalMatrix[0][0] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[1][0] = 0;
    SurrogateNormalMatrix[0][1] = 0;
    SurrogateNormalMatrix[1][1] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[0][2] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[1][2] = fe.surface()->normal().data()[0];
#endif
#if (DIM == 3)
    SurrogateNormalMatrix[0][0] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[1][0] = 0;
    SurrogateNormalMatrix[2][0] = 0;
    SurrogateNormalMatrix[0][1] = 0;
    SurrogateNormalMatrix[1][1] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[2][1] = 0;
    SurrogateNormalMatrix[0][2] = 0;
    SurrogateNormalMatrix[1][2] = 0;
    SurrogateNormalMatrix[2][2] = fe.surface()->normal().data()[2];

    SurrogateNormalMatrix[0][3] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[1][3] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[2][3] = 0;
    SurrogateNormalMatrix[0][4] = fe.surface()->normal().data()[2];
    SurrogateNormalMatrix[1][4] = 0;
    SurrogateNormalMatrix[2][4] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[0][5] = 0;
    SurrogateNormalMatrix[1][5] = fe.surface()->normal().data()[2];
    SurrogateNormalMatrix[2][5] = fe.surface()->normal().data()[1];

#endif
  }

  void CalcTrueNormalMatrix(const ZEROPTV TrueNormal, double (&TrueNormalMatrix)[DIM][3 * (DIM - 1)])
  {
#if (DIM == 2)
    TrueNormalMatrix[0][0] = TrueNormal.x();
    TrueNormalMatrix[1][0] = 0;
    TrueNormalMatrix[0][1] = 0;
    TrueNormalMatrix[1][1] = TrueNormal.y();
    TrueNormalMatrix[0][2] = TrueNormal.y();
    TrueNormalMatrix[1][2] = TrueNormal.x();
#endif
#if (DIM == 3)
    TrueNormalMatrix[0][0] = TrueNormal(0);
    TrueNormalMatrix[1][0] = 0;
    TrueNormalMatrix[2][0] = 0;
    TrueNormalMatrix[0][1] = 0;
    TrueNormalMatrix[1][1] = TrueNormal(1);
    TrueNormalMatrix[2][1] = 0;
    TrueNormalMatrix[0][2] = 0;
    TrueNormalMatrix[1][2] = 0;
    TrueNormalMatrix[2][2] = TrueNormal(2);

    TrueNormalMatrix[0][3] = TrueNormal(1);
    TrueNormalMatrix[1][3] = TrueNormal(0);
    TrueNormalMatrix[2][3] = 0;
    TrueNormalMatrix[0][4] = TrueNormal(2);
    TrueNormalMatrix[1][4] = 0;
    TrueNormalMatrix[2][4] = TrueNormal(0);
    TrueNormalMatrix[0][5] = 0;
    TrueNormalMatrix[1][5] = TrueNormal(2);
    TrueNormalMatrix[2][5] = TrueNormal(1);

#endif
  }

#pragma mark generating terms for integrations

  void CalcCmatrix(double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)])
  {
    /*
     * 3D do not have plane stress case
     */
#if (DIM == 2)
    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      // C for plane stress
      Cmatrix[0][0] = young / (1 - pow(poisson, 2));
      Cmatrix[0][1] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[1][1] = young / (1 - pow(poisson, 2));
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (2 * (1 + poisson));
    }
#endif

    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      // C for plane strain
#if (DIM == 2)
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
#endif
#if (DIM == 3)
      /*
       * this formulation is from FEM book page 241
       */
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      // [previous bug here]
      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      // C for lame parameters
#if (DIM == 2)
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
#endif
#if (DIM == 3)
      double young = mu * (3 * lamda + 2 * mu) / (mu + lamda);
      double poisson = lamda / 2 / (lamda + mu);

      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      // [previous bug here]
      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }
  }

  void CalcBasisMatrix(const TALYFEMLIB::FEMElm &fe, std::vector<std::vector<double>> &Basis)
  {
    assert(Basis.size() == DIM * fe.nbf());
    for (int j = 0; j < fe.nbf(); j++)
    {
#if (DIM == 2)
      /// x-dir
      Basis[2 * j].resize(2);
      /// y-dir
      Basis[2 * j + 1].resize(2);

      Basis[2 * j][0] = fe.N(j);
      Basis[2 * j + 1][0] = 0;
      Basis[2 * j][1] = 0;
      Basis[2 * j + 1][1] = fe.N(j);

#endif
#if (DIM == 3)
    /// x-dir
    Basis[3 * j].resize(3);
    /// y-dir
    Basis[3 * j + 1].resize(3);
    /// z-dir
    Basis[3 * j + 2].resize(3);

    Basis[3 * j][0] = fe.N(j);
    Basis[3 * j + 1][0] = 0;
    Basis[3 * j + 2][0] = 0;

    Basis[3 * j][1] = 0;
    Basis[3 * j + 1][1] = fe.N(j);
    Basis[3 * j + 2][1] = 0;

    Basis[3 * j][2] = 0;
    Basis[3 * j + 1][2] = 0;
    Basis[3 * j + 2][2] = fe.N(j);
#endif
    }
    }

  void CalcBe(const TALYFEMLIB::FEMElm &fe, std::vector<std::vector<double>> &Be)
  {
    assert(Be.size() == DIM * fe.nbf());
    for (int j = 0; j < fe.nbf(); j++)
    {
#if (DIM == 2)
      /// x-dir
      Be[2 * j].resize(3);
      /// y-dir
      Be[2 * j + 1].resize(3);

      Be[2 * j][0] = fe.dN(j, 0);
      Be[2 * j + 1][0] = 0;
      Be[2 * j][1] = 0;
      Be[2 * j + 1][1] = fe.dN(j, 1);
      Be[2 * j][2] = fe.dN(j, 1);
      Be[2 * j + 1][2] = fe.dN(j, 0);
#endif
#if (DIM == 3)
      /// x-dir
      Be[3 * j].resize(6);
      /// y-dir
      Be[3 * j + 1].resize(6);
      /// z-dir
      Be[3 * j + 2].resize(6);

      Be[3 * j][0] = fe.dN(j, 0);
      Be[3 * j + 1][0] = 0;
      Be[3 * j + 2][0] = 0;

      Be[3 * j][1] = 0;
      Be[3 * j + 1][1] = fe.dN(j, 1);
      Be[3 * j + 2][1] = 0;

      Be[3 * j][2] = 0;
      Be[3 * j + 1][2] = 0;
      Be[3 * j + 2][2] = fe.dN(j, 2);

      Be[3 * j][3] = fe.dN(j, 1);
      Be[3 * j + 1][3] = fe.dN(j, 0);
      Be[3 * j + 2][3] = 0;

      Be[3 * j][4] = fe.dN(j, 2);
      Be[3 * j + 1][4] = 0;
      Be[3 * j + 2][4] = fe.dN(j, 0);

      Be[3 * j][5] = 0;
      Be[3 * j + 1][5] = fe.dN(j, 2);
      Be[3 * j + 2][5] = fe.dN(j, 1);
#endif
    }
  }

  void CalcBeCmatrix(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &Be, const double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)], std::vector<std::vector<double>> &BeCmatrix)
  {
    assert(BeCmatrix.size() == DIM * fe.nbf());
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      BeCmatrix[a].resize(3 * (DIM - 1));
      for (int b = 0; b < 3 * (DIM - 1); b++)
      {
        double N = 0;
        for (int k = 0; k < 3 * (DIM - 1); k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        BeCmatrix[a][b] = N;
      }
    }
  }

  void CalcStressDotNormal(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &BeCmatrix, const double (&NormalMatrix)[DIM][3 * (DIM - 1)], std::vector<std::vector<double>> &StressDotNormal)
  {
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      StressDotNormal[a].resize(DIM);
      for (int b = 0; b < DIM; b++)
      {
        for (int k = 0; k < 3 * (DIM - 1); k++)
        {
          StressDotNormal[a][b] += BeCmatrix[a][k] * NormalMatrix[b][k];
        }
      }
    }
  }

  void CalcBeCmatrix(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &Be, const std::vector<std::vector<double>> Cmatrix, std::vector<std::vector<double>> &BeCmatrix)
  {
    assert(BeCmatrix.size() == DIM * fe.nbf());
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      // BeCmatrix[a].resize(3*(DIM-1));
      for (int b = 0; b < 3 * (DIM - 1); b++)
      {
        double N = 0;
        for (int k = 0; k < 3 * (DIM - 1); k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        BeCmatrix[a][b] = N;
      }
    }
  }
  void CalcGradStrainDotDistance(const TALYFEMLIB::FEMElm &fe, std::vector<std::vector<double>> &GradStressDotDistance, const double (&distance)[DIM])
    {
        // the shape of this matrix is similar as the BeCmatrix
//        this is the number of rows
        assert(GradStressDotDistance.size() == DIM * fe.nbf());
//        for each row have columns = 3*(DIM-1)

#if (DIM==2)
        for(int a=0; a<fe.nbf(); a++)
    {
        GradStressDotDistance[2*a].resize(3);  // 3 is the number of columns in BMatrix
        GradStressDotDistance[2*a+1].resize(3); // 3 is the number of columns in BMatrix
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y : for first column
        GradStressDotDistance[2*a][0] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1];
        GradStressDotDistance[2*a][1] = 0;
        GradStressDotDistance[2*a][2] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1];
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y : for second column
        GradStressDotDistance[2*a+1][0] = 0;
        GradStressDotDistance[2*a+1][1] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1];
        GradStressDotDistance[2*a+1][2] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1];

    }
#endif
#if (DIM==3)
        for(int a=0; a<fe.nbf(); a++)
    {
//        print the value of the second derivative of the shape function



        GradStressDotDistance[3*a].resize(6);  // 6 is the number of columns in BMatrix
        GradStressDotDistance[3*a+1].resize(6); // 6 is the number of columns in BMatrix
        GradStressDotDistance[3*a+2].resize(6); // 6 is the number of columns in BMatrix
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for first column
        GradStressDotDistance[3*a][0] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a][1] = 0;
        GradStressDotDistance[3*a][2] = 0;
        GradStressDotDistance[3*a][3] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1] + fe.d2N(a,1,2)*distance[2];
        GradStressDotDistance[3*a][4] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
        GradStressDotDistance[3*a][5] = 0;
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for second column
        GradStressDotDistance[3*a+1][0] = 0;
        GradStressDotDistance[3*a+1][1] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1] + fe.d2N(a,1,2)*distance[2];
        GradStressDotDistance[3*a+1][2] = 0;
        GradStressDotDistance[3*a+1][3] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a+1][4] = 0;
        GradStressDotDistance[3*a+1][5] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for third column
        GradStressDotDistance[3*a+2][0] = 0;
        GradStressDotDistance[3*a+2][1] = 0;
        GradStressDotDistance[3*a+2][2] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
        GradStressDotDistance[3*a+2][3] = 0;
        GradStressDotDistance[3*a+2][4] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a+2][5] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];

    }
#endif
// Check  if any of the value is nan

    }
  DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe)
  {
    return pow((pow(2, DIM) * fe.volume_jacc()), (double)1 / DIM);
  }
  ZEROPTV CalcNormalTraction(const TALYFEMLIB::FEMElm &fe, const double d[DIM])
  {
//    COMPUTE THE TRUE POSITION
    double eps = 1e-10;
    ZEROPTV traction = ZEROPTV{0, 0, 0};
    idata_->traction_vector_;


    if (fabs(fe.position().x()-idata_->mesh_def.physDomain.max[0]) < eps)
    {

#ifndef NDEBUG
//        CREATE A FILE TO WRITE THE TRACTION
        std::ofstream file;
        file.open("traction.csv", std::ios::app);
//        IF THE FILE IS EMPTY THEN WRITE THE HEADER
        if (file.tellp() == 0)
        {
          file << "X,Y,Z" << std::endl;
        }

//        this is appending right ???
        file << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << ",\n"<< std::endl;
        file.close();

#endif
//print the traction
      traction = ZEROPTV{idata_->NormalTraction.traction_vec[0],idata_->NormalTraction.traction_vec[1], idata_->NormalTraction.traction_vec[2]};
      return traction;
    }
    if(fabs(fe.position().y()-idata_->mesh_def.physDomain.max[1]) < eps)
    {
      traction = ZEROPTV{0, 0, 0};
      return traction;
    }
    return traction;

  }
  ZEROPTV CalcNormalTractionUVL(const TALYFEMLIB::FEMElm &fe, const double h)
    {
        double eps = 1e-10;
        ZEROPTV traction = ZEROPTV{0, 0, 0};
        if (fabs(fe.position().x() - idata_->mesh_def.physDomain.max[0]) < eps)
        {
            double y = fe.position().y();
            double max_y = idata_->mesh_def.physDomain.max[1];
            double min_y = idata_->mesh_def.physDomain.min[1];

//            at this point applying the traction shifts the BC to the left
//            this is not the part of the domain and doesn't require additional traction
            if(y >= 1.0)
            {
                return traction;
            }

            traction = ZEROPTV{1, 0, 0};
            return traction;

//            double y = fe.position().y();
//            double max_y = idata_->mesh_def.physDomain.max[1];
//            double min_y = idata_->mesh_def.physDomain.min[1];
//
//// Compute the mid-point of the domain along the y-axis
//            double mid_y = (max_y + min_y) / 2.0;
//
//// Compute the parabolic profile based on the y-position
//            double parabolic_profile = 1.0 - pow((y - mid_y) / (max_y - mid_y), 2);
//
//// Apply the traction with a parabolic distribution
//            traction = ZEROPTV{parabolic_profile, 0, 0};
//
//            std::cout << "Traction for Query Point:" << std::endl;
//            std::cout << traction << std::endl;
//            return traction;
        }
        if (fabs(fe.position().y() - idata_->mesh_def.physDomain.max[1]) < eps)
        {
            traction = ZEROPTV{0, 0, 0};
            return traction;
        }
        return traction;
    }
  ZEROPTV computeTraction(const TALYFEMLIB::FEMElm& fe, const LEInputData* idata_) {

//    cantilever beam normal force:

//
        assert(idata_ != nullptr);  // Make sure input data is not null

        size_t num_results = 2;  // Searching for 2 nearest neighbors
        std::vector<uint32_t> ret_index(num_results);
        std::vector<double> out_dist_sqr(num_results);

#if (DIM == 2)
        const double query_pt[3] = {fe.position().x() - idata_->CsvForce.shift_pos[0], fe.position().y() - idata_->CsvForce.shift_pos[1], 0.0};
#endif

#if (DIM == 3)
        const double query_pt[3] = {fe.position().x() - idata_->CsvForce.shift_pos[0], fe.position().y() - idata_->CsvForce.shift_pos[1], fe.position().z() - idata_->CsvForce.shift_pos[2]};
#endif

        num_results = idata_->kdTree_->knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        auto closest_point_1 = idata_->traction_position_.pts[ret_index[0]];
        auto closest_point_2 = idata_->traction_position_.pts[ret_index[1]];

        ZEROPTV traction1 = idata_->traction_vector_[ret_index[0]];
        ZEROPTV traction2 = idata_->traction_vector_[ret_index[1]];

        ZEROPTV traction;

#if (DIM == 2)
        traction = util_funcs::linear_interpolation_extrapolation(ZEROPTV{closest_point_1.x, closest_point_1.y, 0.0}, ZEROPTV{closest_point_2.x, closest_point_2.y, 0.0}, traction1, traction2, ZEROPTV{query_pt[0], query_pt[1], 0.0});
#endif

#if (DIM == 3)
        traction = util_funcs::linear_interpolation_extrapolation(ZEROPTV{closest_point_1.x, closest_point_1.y, closest_point_1.z}, ZEROPTV{closest_point_2.x, closest_point_2.y, closest_point_2.z}, traction1, traction2, ZEROPTV{query_pt[0], query_pt[1], query_pt[2]});
#endif

#ifndef NDEBUG

#endif
              std::cout << "Query Point:" << std::endl;
              std::cout << "X: " << query_pt[0] << std::endl;
              std::cout << "Y: " << query_pt[1] << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << query_pt[2] << std::endl;
#endif
              std::cout << "Traction for Query Point:" << std::endl;
              std::cout << traction << std::endl;


              std::cout << "1st Closest Point:" << std::endl;
              std::cout << "X: " << closest_point_1.x << std::endl;
              std::cout << "Y: " << closest_point_1.y << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << closest_point_1.z << std::endl;
#endif
              std::cout << "Traction1 for 1st Closest Point:" << std::endl;
              std::cout << traction1 << std::endl;


              std::cout << "2nd Closest Point:" << std::endl;
              std::cout << "X: " << closest_point_2.x << std::endl;
              std::cout << "Y: " << closest_point_2.y << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << closest_point_2.z << std::endl;
#endif
              std::cout << "Traction2 for 2nd Closest Point:" << std::endl;
              std::cout << traction2 << std::endl;

        return traction;
    }



//    some utility function based on requirements
    bool GpOnDomainWall(const TALYFEMLIB::FEMElm &fe)
    {
        /// All the possible cases
        static const double eps = 1e-7;

        bool x_minus_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.min[0]) < eps;
        bool y_minus_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.min[1]) < eps;
        bool x_max_wall = fabs(fe.position().x() - idata_->mesh_def.physDomain.max[0]) < eps;
        bool y_max_wall = fabs(fe.position().y() - idata_->mesh_def.physDomain.max[1]) < eps;
        bool z_minus_wall = false;
        bool z_max_wall = false;
#if(DIM == 3)
        z_minus_wall = fabs(fe.position().z() - idata_->mesh_def.physDomain.min[2]) < eps;
        z_max_wall = fabs(fe.position().z() - idata_->mesh_def.physDomain.max[2]) < eps;
#endif

        return x_minus_wall || y_minus_wall || x_max_wall || y_max_wall || z_minus_wall || z_max_wall;
    }

    void saveCoordinatesToFile(const TALYFEMLIB::FEMElm &fe, const string &filename) {
        // Open the file in append mode
        std::ofstream file(filename, std::ios::app);

        // Check if the file is open
        if (file.is_open()) {
            // Write the coordinates to the file
            file << fe.position().x() << ", " << fe.position().y() << ", " << fe.position().z() << "\n";

            // Close the file
            file.close();
        } else {
            std::cerr << "Unable to open file for writing\n";
        }
    }

};

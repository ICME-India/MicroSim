    //! Elastic stress, strain, displacement fields only for precipitate growth
    iCorr=0;
    
    if (swch == 2)
    {    
    if (phases == 3)
    {
            do
            {
                 
            {
                
                fvVectorMatrix DEqn
                (
                    dimt*dimt*fvm::d2dt2(D)
                 ==
                    dimx*dimx*fvm::laplacian(2*(mu2_elast*phi_3*phi_3*(3-2*phi_3) + mu1_elast*(1-phi_3)*(1-phi_3)*(1+2*phi_3)) 
		    + lambda2*phi_3*phi_3*(3-2*phi_3) + lambda1*(1-phi_3)*(1-phi_3)*(1+2*phi_3), D, "laplacian(DD,D)")
                  + dimx*dimx*divSigmaExp
		  -dimx*dimx*fvc::div((2*mu1_elast*(1-phi_3)*(1-phi_3)*(1+2*phi_3) + 2*mu2_elast*phi_3*phi_3*(3-2*phi_3))*(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1
                             + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)
                             +(lambda2*phi_3*phi_3*(3-2*phi_3) + (1-phi_3)*(1-phi_3)*(1+2*phi_3)*lambda1)*I*tr(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1
                             + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2))
                );

                /*
                     DEqn -= (sig3)*fvc::div((2*mu1*(1-phi_4)*(1-phi_4)*(1+2*phi_4) + 2*mu2*phi_4*phi_4*(3-2*phi_4))*(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 
			     + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)
			     +(lambda2*phi_4*phi_4*(3-2*phi_4) + (1-phi_4)*(1-phi_4)*(1+2*phi_4)*lambda1)*I*tr(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 
			     + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3));
		*/

                InitialResidual_3 = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                
                D.correctBoundaryConditions();
                gradD = fvc::grad(D);
		//gradD = gradD.T();
    /*
		scalar eigen_xx=0.0; 
		scalar eigen_yy=0.0;
		scalar eigen_zz=0.0;
		    
                  forAll(strain,cellI)
                  {   
		    eigen_xx  =   phi_1[cellI]*phi_1[cellI]*(3-2*phi_1[cellI])*cEigenStrain1.component(symmTensor::XX).value();
		    eigen_xx +=   phi_2[cellI]*phi_2[cellI]*(3-2*phi_2[cellI])*cEigenStrain2.component(symmTensor::XX).value();
		    eigen_xx +=   phi_3[cellI]*phi_3[cellI]*(3-2*phi_3[cellI])*cEigenStrain3.component(symmTensor::XX).value();
		    
		    eigen_yy  =   phi_1[cellI]*phi_1[cellI]*(3-2*phi_1[cellI])*cEigenStrain1.component(symmTensor::YY).value();
		    eigen_yy +=   phi_2[cellI]*phi_2[cellI]*(3-2*phi_2[cellI])*cEigenStrain2.component(symmTensor::YY).value();
		    eigen_yy +=   phi_3[cellI]*phi_3[cellI]*(3-2*phi_3[cellI])*cEigenStrain3.component(symmTensor::YY).value();
		    
		    eigen_zz  =   phi_1[cellI]*phi_1[cellI]*(3-2*phi_1[cellI])*cEigenStrain1.component(symmTensor::ZZ).value();
		    eigen_zz +=   phi_2[cellI]*phi_2[cellI]*(3-2*phi_2[cellI])*cEigenStrain2.component(symmTensor::ZZ).value();
		    eigen_zz +=   phi_3[cellI]*phi_3[cellI]*(3-2*phi_3[cellI])*cEigenStrain3.component(symmTensor::ZZ).value();
		    
                    strain[cellI].component(symmTensor::XX) = gradD[cellI].component(tensor::XX) - eigen_xx;
                    strain[cellI].component(symmTensor::YY) = gradD[cellI].component(tensor::YY) - eigen_yy;
                    strain[cellI].component(symmTensor::ZZ) = gradD[cellI].component(tensor::ZZ) - eigen_zz;
                    strain[cellI].component(symmTensor::XY) = 0;
                    strain[cellI].component(symmTensor::XZ) = 0;
                    strain[cellI].component(symmTensor::YZ) = 0;
                  }*/
		
   strain=((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)))
                &&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
          +((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)))
                &&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
          +((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)))
                &&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);
		
	
                sigmaD = (mu2_elast*phi_3*phi_3*(3-2*phi_3) + mu1_elast*(1-phi_3)*(1-phi_3)*(1+2*phi_3))*twoSymm(gradD) 
		+ (lambda2*phi_3*phi_3*(3-2*phi_3)    + lambda1*(1-phi_3)*(1-phi_3)*(1+2*phi_3))*(I*tr(gradD))
                + (mu2_elast_*phi_3*phi_3*(3-2*phi_3) + mu1_elast_*(1-phi_3)*(1-phi_3)*(1+2*phi_3))*strain;
            
                if (compactNormalStress)
                {
                    divSigmaExp = fvc::div
                    (
                        sigmaD - (2*mu2_elast*phi_3*phi_3*(3-2*phi_3) + 2*mu1_elast*(1-phi_3)*(1-phi_3)*(1+2*phi_3)
			+lambda2*phi_3*phi_3*(3-2*phi_3) + (1-phi_3)*(1-phi_3)*(1+2*phi_3)*lambda1)*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += fvc::div(sigmaD);
                }
            }

        } while (InitialResidual_3 > convergenceTolerance && ++iCorr < nCorr);

    	Info<< "Min/max Dx:" << min(D.component(0)).value() << ' ' << max(D.component(0)).value() << endl;
	Info<< "Min/max Dy:" << min(D.component(1)).value() << ' ' << max(D.component(1)).value() << endl;
	if (dimensions == 3)
	{
	Info<< "Min/max Dz:" << min(D.component(2)).value() << ' ' << max(D.component(2)).value() << endl;
	}
    }
    
    }

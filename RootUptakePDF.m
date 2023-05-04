function Fr=RootUptakePDF(Theta,z,Zr,type,Seg_depth)
     nvar=length(Theta);
% %--------------------------------------------------------------------------
%                  one parameter root water uptake models
%--------------------------------------------------------------------------
aa=Theta(1);  
%  if nvar==2                      % one parameter root water uptake model

 if nvar==3
     bb=Theta(2); end
 if nvar==4
     bb=Theta(2); cc=Theta(3);
 end
         switch type
            case '3MixUniform'
                d1=Seg_depth(1);                % First segment
                d2=Seg_depth(2)+d1;             % Second segment
               Fr= aa*unifpdf(z,0,d1)+bb*unifpdf(z,d1,d2)+(1-aa-bb)*unifpdf(z,d2,Zr);
%                Fr=aa*(z >= 0 && z < d1)+bb*(z>=d1 && z<d2)+(1-aa*d1-bb*(d2-d1))*(z>=d2 &&z<=Zr);
            case 'UnitLindley'                  % 0<x<1  aa>0
                if(z<0 ||z>=Zr)
                    Fr=0;
                else
                    x =z/Zr;
                    Fr=aa^2/(1+aa)/(1-x)^3*exp(-aa*x/(1-x))/Zr;
                end
            case 'Lindley'                      % z>0    aa>0
                if(z<=0)
                    Fr=0;
                else
                    Fr=aa^2/(1+aa)*(1+z)*exp(-aa*z);
                end
            case 'TrExponent'       %Truncated Exponential distribution
                if(z<0 ||z >Zr)
                    Fr=0;
                else
                    y=z/Zr;
                    if aa== 0
                        Fr=1/Zr;
                    else
                        Fr=aa*exp(-aa*y)/(1-exp(-aa))/Zr;  %Cook model
                    end
                end
            case 'Pert1'            % one parameter Pert model
                aa1=1+4*(aa/Zr); bb1=1+4*(Zr-aa)/Zr;
                lnFr=(aa1-1)*log(z)+(bb1-1)*log(Zr-z)-betaln(aa1,bb1)-(aa1+bb1-1)*log(Zr);
                Fr=exp(lnFr);
        
%--------------------------------------------------------------------------
%                  Two one-parameter root water uptake models
%--------------------------------------------------------------------------
      
            case 'MixLindley'
                if(z<=0)
                    Fr=0;
                else
                    Fr1=@(aa1) aa1^2/(1+aa1)*(1+z)*exp(-aa1*z);
                    Fr=cc*Fr1(aa)+(1-cc)*Fr1(bb);
                end
            case 'MixULindley'
                if(z<=0 ||z>=Zr)
                    Fr=0;
                else
                    x=z/Zr;
                    Fr1=aa^2/(1+aa)/(1-x)^3*exp(-aa*x/(1-x))/Zr;
                    Fr2=bb^2/(1+bb)/(1-x)^3*exp(-bb*x/(1-x))/Zr;
                    Fr=cc*Fr1+(1-cc)*Fr2;
                end
%             case '3TrWeibull'
%                 if z>cc
%                     Fr=bb/aa*((z-cc)/aa)^(bb-1)*exp(-((z-cc)/z)^bb);
%                 else
%                     Fr=0;
%                 end
            case '3UnitWeibull' % doi.org/10.1080/02664763.2019.1657813
                y=z/Zr-cc;
                if(y<=0 )
                    Fr=0;
                else
                    %Fr=beta/alpha*(z/alpha)^(beta-1)*exp(-(z/alpha)^beta);
                    lnFr=log(bb)+log(aa)+(bb-1)*log(-log(y))-aa*(-log(y))^bb;
                    Fr=exp(lnFr)/y/Zr;
                end
            case '3Lognormal'
                y=z/Zr-cc;
                if(y<=0 )
                    Fr=0;
                else
                    Phi= 1/y/bb/sqrt(2*3.14)*exp(-0.5*(log(y/aa)/bb)^2);
                    BigPhi= 0.5*(1+erf(log((1-cc)/aa)/bb/sqrt(2)));
                    Fr=Phi/BigPhi/Zr;
                end
            case '3Tr2SNormal'
                if z>=Zr|| z<=0
                    Fr=0;
                else
                    y=@(x)2/cc*sinh(x/bb-aa/bb);
                    Phi= 1/sqrt(2*3.14)*exp(-0.5*y(z)*y(z));
                    BigPhi=@ (x) 0.5*(1+erf(x/sqrt(2)));
                    Fr=2/cc/bb*cosh((z-aa)/bb)*Phi;
                    cdf_D=BigPhi(y(Zr))-BigPhi(y(0));
                    Fr=Fr/cdf_D;
                end
            case '3Bimodel'
                y=(z-aa)/bb;
                Phi= 1/bb/sqrt(6.28)*exp(-0.5*y^2);
                BigPhi=@ (x) 0.5*(1+erf(x/sqrt(2)))+cc*(2-cc*x)/(2+cc^2)*1/sqrt(6.28)*exp(-0.5*x^2);
                Fr=Phi*((1+cc*y)^2+1)/(2+cc^2)/(BigPhi((Zr-aa)/bb)-BigPhi(-aa/bb));

       
%--------------------------------------------------------------------------
%                  Two-parameter root water uptake models
%--------------------------------------------------------------------------
          
            case 'QuasiTwoPaLindley'  %z>0 aa>-1. When bb=aa, it is Lindley. When bb=0, it is gamma
                Fr=aa*(bb+aa*z)/(bb+1)*exp(-aa*z);
            
            case 'UnitInvGaus'   %bb/2/aa/aa>=1;
                if z<=0 || z>=Zr
                    Fr=0;
                else
                    x=z/Zr;
                    Fr=sqrt(bb/2/3.14)/x/(-log(x))^1.5*exp(bb/2/aa/aa/log(x)*(log(x)+aa)^2)/Zr;
                end
            case 'ToppLeone'
                x=z/Zr;
                Fr=2*aa*bb*(1-x)*(x*(2-x))^(aa-1)*(1-(x*(2-x))^aa)^(bb-1)/Zr;
            case 'JohnsonSb'
                if z<=0 || z>=Zr
                    Fr=0;
                else
                    x=z/Zr;
                    Fr=bb/sqrt(2*3.14)/x/(1-x)*exp(-0.5*(aa+bb*log(x/(1-x)))^2)/Zr;
                end
            case 'Simplex'
                if z<=0 || z>=Zr
                    Fr=0;
                else
                    x=z/Zr;
                    dx=(x-aa)^2/x/(1-x)/(1-aa)^2/aa/aa;
                    Fr=1/sqrt(2*3.14*bb*(x*(1-x))^3)*exp(-dx/2/bb)/Zr;
                end
            case 'UnitGompertz'
                if z<=0 || z>=Zr
                    Fr=0;
                else
                    x=z/Zr;
                    Fr=aa*bb/x^(bb+1)*exp(-aa/x^bb+aa)/Zr;
                end
            case 'TrWeibull'
                if(z<=0 )
                    Fr=0;
                else
                    %Fr=beta/alpha*(z/alpha)^(beta-1)*exp(-(z/alpha)^beta);
                    lnFr=log(bb)-log(aa)+(bb-1)*log(z/aa)-(z/aa)^bb;
                    Weibcdf=1-exp(-(Zr/aa)^bb);
                    Fr=exp(lnFr)/Weibcdf;
                end
            case 'UnitWeibull' % doi.org/10.1080/02664763.2019.1657813
                y=z/Zr;
                if(z<=0 )
                    Fr=0;
                else
                    %Fr=beta/alpha*(z/alpha)^(beta-1)*exp(-(z/alpha)^beta);
                    lnFr=log(bb)+log(aa)+(bb-1)*log(-log(y))-aa*(-log(y))^bb;
                    Fr=exp(lnFr)/y/Zr;
                end

            case 'Beta'
                if(z<=0 ||z>=Zr)
                    Fr=0;
                else
                    lnFr = (aa-1)*log(z/Zr)+(bb-1)*log(1-z/Zr)-betaln(aa,bb)-log(Zr);
                    Fr = exp(lnFr);
                end
            case 'TrGamma'
                kap=aa*aa/bb/bb; lmd=bb*bb/aa;
                if(z<=0 )
                    Fr=0;
                else
                    lnFr=-gammaln(kap)-kap*log(lmd)+(kap-1)*log(z)-z/lmd;
                    Fr=exp(lnFr)/gamcdf(Zr,kap,lmd);
                end
            case 'UnitGamma1'   %doi.org/10.1111/j.1467-842X.1977.tb01277.x
                kap=aa*aa/bb/bb; lmd=bb*bb/aa;
                x=z/Zr;
                if(z<=0||z>=Zr )
                    Fr=0;
                else
                    lnFr=-gammaln(kap)-kap*log(lmd)+(kap-1)*log(-log(x))+log(x)*(1/lmd-1);
                    %Fr=aa^bb/gamma(bb)*(z/Zr)^(bb-1)*(log(Zr/z))^(bb-1)/Zr;
                    Fr=exp(lnFr)/Zr;
                end
            case 'UnitGamma2'
                kap=aa*aa/bb/bb; lmd=bb*bb/aa;
                x=z/Zr;
                if(z<=0||z>=Zr )
                    Fr=0;
                else
                    lnFr=-gammaln(kap)-kap*log(lmd)+(kap-1)*log(-log(1-x))+log(1-x)*(1/lmd);
                    Fr=exp(lnFr)/x/Zr;
                end
            case 'Kumaraswamy'
                if(z<=0 ||z>=Zr)
                    Fr=0;
                else
                    %Fr=alpha*beta*(z/Zr)^(alpha-1)*(1-(z/Zr)^alpha)^(beta-1)/Zr;
                    lnFr=log(aa)+log(bb)+(aa-1)*log(z/Zr)+(bb-1)*log(1-(z/Zr)^aa)-log(Zr);
                    Fr=exp(lnFr);
                end
                %------------------------------------------------------------------------
                %The Cauchy distribution has a very heavy tail, comparable to the tail of
                % the Pareto (1, c) distribution.   In fact, the tail is so heavy that the
                % distribution does not have a mean value. ∫fdx=1 but ∫xfdx does not exist
                % and so  % the mean of X does not exist.
                %------------------------------------------------------------------------
            case 'TrCauchy'
                Fr=1/3.14/bb/(1+((z-aa)/bb)^2);
                CDF_diff=1/3.14*(atan((Zr-aa)/bb)-atan(-aa/bb));
                Fr=Fr/CDF_diff;
                %-------------------------------------------------------------------------
                %The PERT distribution also uses the most likely value, but it is designed
                % to generate a distribution that more closely resembles realistic
                % probability distribution. Depending on the values provided, the PERT
                % distribution can provide a close fit to the normal or lognormal
                % distributions.Like the triangular distribution, the PERT distribution
                % emphasizes the “most likely” value over the minimum and maximum estimates.
                % However, unlike the triangular distribution the PERT distribution
                % constructs a smooth curve which places progressively more emphasis on
                % values around (near) the most likely value, in favor of values around
                % the edges. In practice, this means that we “trust” the estimate for
                % the most likely value, and we believe that even if it is not exactly
                % accurate (as estimates seldom are), we have an expectation that the
                % resulting value will be close to that estimate.

            case 'Pert2'
                aa1=1+bb*(aa/Zr); bb1=1+bb*(Zr-aa)/Zr;
                lnFr=(aa1-1)*log(z)+(bb1-1)*log(Zr-z)-betaln(aa1,bb1)-(aa1+bb1-1)*log(Zr);
                Fr=exp(lnFr);
                %--------------------------------------------------------------------------
                % It is a symmetrical distribution, unimodal (it has one peak) and is similar
                % in shape to the normal distribution. In fact, the logistic and normal
                % distributions are so close in shape (although the logistic tends to
                % have slightly fatter tails) that for most applications it’s impossible
                % to tell one from the other. The logistic distribution is mainly used
                % because the curve has a relatively simple cumulative distribution
                % formula to work with. The formula approximates the normal distribution
                % extremely well.
                %
            case 'TrLogistic'
                Fr=0.25/bb*(sech((z-aa)/2/bb))^2;
                CDF_diff=0.5*tanh((Zr-aa)/2/bb)-0.5*tanh(-aa/2/bb);
                Fr=Fr/CDF_diff;
                %-------------------------------------------------------------------------
                %It is similar in shape to the log-normal distribution but has heavier tails.
                % Unlike the log-normal, its cumulative distribution function can be
                % written in closed form. The parameter α > 0 is a scale parameter
                % and is also the median of the distribution. The parameter β > 0
                % is a shape parameter. The distribution is unimodal when β > 1
                % and its dispersion decreases as β  increases.
            case 'LogLogistic'
                if z<=0
                    Fr=0;
                else
                    CDF=1/(1+(Zr/aa)^(-bb));
                    Fr=bb/aa*(z/aa)^(bb-1)/(1+(z/aa)^bb)^2;
                    Fr=Fr/CDF;
                end
            case 'UnitLogistic'  % aa is the median. doi: 10.1590/0101-7438.2018.038.03.0555
                if (z<=0 ||z>=Zr)
                    Fr=0;
                else
                    x=z/Zr;  aa=aa/Zr;
                    LnFr=log(bb)+bb*log(aa)+(bb-1)*log(x*(1-x))+bb*log(1-aa)-...
                        2*log(x^bb*(1-aa)^bb+aa^bb*(1-x)^bb);
                    Fr=exp(LnFr)/Zr;
                end
            case 'Lognormal'
                if(z<=0 )
                    Fr=0;
                else
                    Phi= 1/z/bb/sqrt(2*3.14)*exp(-0.5*(log(z/aa)/bb)^2);
                    BigPhi= 0.5*(1+erf(log(Zr/aa)/bb/sqrt(2)));
                    Fr=Phi/BigPhi;
                end

            case 'Tr2Normal'
                Phi= 1/bb/sqrt(2*3.14)*exp(-0.5*((z-aa)/bb)^2);
                BigPhi=@ (x) 0.5*(1+erf(x/sqrt(2)));
                CDF0=BigPhi(-aa/bb); CDFZr=BigPhi((Zr-aa)/bb);
                if(CDF0<0.01 && CDFZr>0.99)
                    Fr=Phi;
                else
                    Fr=Phi/(CDFZr-CDF0);
                end
                %--------------------------------------------------------------------------
                %Recall that the logit is the log odds of a probability. So for starters
                % we know that the logit-normal has the domain as a probability, being
                % bound between 0 and 1 (technically a probability can be 0 and 1 exactly,
                % where as the logit-normal is not defined for these extremes). The odds
                % are the ratio of the probability of one event happening over an other,
                % in the case of the logit, we are comparing the probability of an event
                % happening, pp to the probability of it not happening, 1−p.
                % o=p1/(1−p​)​​
                % For the log part of logit we just perform a log transform of our odds.
                % So we can mathematically define logit as: logit(p)=log(p/(1-p)).
            case 'LogitNorm'
                if z<=0 || z>=Zr
                    Fr=0;
                else
                    zx=(log(z/(Zr-z))-log(aa/(Zr-aa)))/bb;
                    Fr=1/bb/sqrt(6.28)*exp(-0.5*zx^2)*Zr/z/(Zr-z);
                end
                %-------------------------------------------------------------------------
                %Given a normally distributed random variable X with mean μ and variance σ2,
                % the random variable Y = |X| has a folded normal distribution. Such a case
                % may be encountered if only the magnitude of some variable is recorded,
                % but not its sign. The distribution is called "folded" because probability
                % mass to the left of x = 0 is folded over by taking the absolute value.
                % In the physics of heat conduction, the folded normal distribution is a
                % fundamental solution of the heat equation on the half space; it
                % corresponds to having a perfect insulator on a hyperplane through the origin.
            case 'FoldedNorm'
                y1=(z-aa)/bb;           y2=(z+aa)/bb;
                Phi1= 1/bb/sqrt(6.28)*exp(-0.5*y1^2);
                Phi2= 1/bb/sqrt(6.28)*exp(-0.5*y2^2);
                CDF=0.5*erf((Zr+aa)/sqrt(2)/bb)+0.5*erf((Zr-aa)/sqrt(2)/bb);
                Fr=(Phi1+Phi2)/CDF;
            case 'UnitFoldN'
                if (z<=0 || z>=Zr)
                    Fr=0;
                else
                    y1=(atanh(z/Zr)-aa)/bb;           y2=(atanh(z/Zr)+aa)/bb;
                    Phi1= 1/bb/sqrt(6.28)*exp(-0.5*y1^2);
                    Phi2= 1/bb/sqrt(6.28)*exp(-0.5*y2^2);
                    Fr=1/(1-(z/Zr)^2)*(Phi1+Phi2)/Zr;
                end
            case 'Tr2SNormal'
                if z>=Zr|| z<=0
                    Fr=0;
                else
                    cc=2;
                    y=@(x)2/cc*sinh(x/bb-aa/bb);
                    Phi= 1/sqrt(2*3.14)*exp(-0.5*y(z)*y(z));
                    BigPhi=@ (x) 0.5*(1+erf(x/sqrt(2)));
                    Fr=2/cc/bb*cosh((z-aa)/bb)*Phi;
                    cdf_D=BigPhi(y(Zr))-BigPhi(y(0));
                    Fr=Fr/cdf_D;
                end
            case 'Bimodel'
                cc=5;
                y=(z-aa)/bb;
                Phi= 1/bb/sqrt(6.28)*exp(-0.5*y^2);
                BigPhi=@ (x) 0.5*(1+erf(x/sqrt(2)))+cc*(2-cc*x)/(2+cc^2)*1/sqrt(6.28)*exp(-0.5*x^2);
                Fr=Phi*((1+cc*y)^2+1)/(2+cc^2)/(BigPhi((Zr-aa)/bb)-BigPhi(-aa/bb));

                %-------------------------------------------------------------------------
                %The Wigner semicircle distribution, named after the physicist Eugene Wigner,
                % is the probability distribution on [−R, R] whose probability density
                % function f is a scaled semicircle (i.e., a semi-ellipse) centered at (0, 0):
                % It is also a scaled beta  % distribution: if Y is beta-distributed with
                % parameters α = β = 3/2, then X = 2RY – R has the Wigner semicircle distribution.
            case 'Wigner'
                if abs(z-aa)>=bb
                    Fr=0;
                else
                    Fr=2/3.14/bb/bb*sqrt(bb*bb-(z-aa).^2);
                end
             case 'CSpline'
                 x=[0, 20, 40 ,80 ,180 ,250];
                 y=[0,aa,bb,cc,0,0];
                 pp=pchip(x,y);
                 CDF=integral(@(xx)ppval(pp,xx),0,Zr);
                 Fr = fnval(pp,z)/CDF;
        end
     end
 

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //setting Defaoult Values
    ui->P1->setText("0.0155"); ui->P2->setText("0.296632"); ui->P3->setText("0.060015"); ui->P4->setText("-0.4515");
    ui->P5->setText("0.296632"); ui->P6->setText("-0.06055"); ui->P7->setText("0.453"); ui->P8->setText("0.0");
    ui->P9->setText("0.001260"); ui->P10->setText("0.0"); ui->P11->setText("7.36");

    ui->R1->setText("0.0015"); ui->R2->setText("0.025"); ui->R3->setText("0.015"); ui->R4->setText("-0.01");
    ui->R5->setText("0.02"); ui->R6->setText("-0.015"); ui->R7->setText("0.075"); ui->R8->setText("0.0");
    ui->R9->setText("0.0"); ui->R10->setText("-0.175"); ui->R11->setText("0.05");

    ui->AOA->setText("5"); ui->uinf->setText("1");



}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::Shape_Function(std::vector<double> a, std::vector<double> b, double X, double &Z_u0, double &Z_d0)
{
    //inputs:Polynomials' coffients , x coordinate
    //output:z_upper coordinate,z_lower coordinate
    double chord=1.0;

    //legendre Polynomials
    double P2=(2*X/chord)-1;
    double P3=(6*pow(X/chord,2)) -(6*(X/chord)) + 1 ;
    double P4=(20*pow(X/chord,3)) - (30*pow(X/chord,2)) + (12*X/chord) -1;
    double P5=(70*pow(X/chord,4)) - (140*pow(X/chord,3)) + (90*pow(X/chord,2)) - (20*X/chord) + 1;
    double P6=(252*pow(X/chord,5)) -(630*pow(X/chord,4)) +(560*pow(X/chord,3)) -(210*pow(X/chord,2)) +(30*X/chord) - 1;

    //the upper and lower z coordinates
    Z_u0= (pow((1-X),3))*(sqrt(a[0]*X)+a[1]*(P2+1)+a[2]*(P3-1)+a[3]*(P4+1)+a[4]*(P5-1)+a[5]*(P6+1));
    Z_d0= -(pow((1-X),3)* (sqrt(b[0]*X)+ b[1]*(P2+1) +(b[2]*(P3-1)) + (b[3]*(P4+1)) +(b[4]*(P5-1)) + (b[5]*(P6+1))));

}

void MainWindow::Shape_Function_PARSIK(std::vector<double> a, double X, double &Z_u0, double &Z_d0)
{
    //inputs:Polynomials' coffients , x coordinate
    //output:z_upper coordinate,z_lower coordinate


    //the upper and lower z coordinates
    Z_u0= a[0]*pow(X,0.5) + a[1]*pow(X,1.5) +a[2]*pow(X,2.5) +a[3]*pow(X,3.5) +a[4]*pow(X,4.5) +a[5]*pow(X,5.5);
    Z_d0= a[6]*pow(X,0.5) + a[7]*pow(X,1.5) +a[8]*pow(X,2.5) +a[9]*pow(X,3.5) +a[10]*pow(X,4.5) +a[11]*pow(X,5.5);



}

void MainWindow::Parsec(std::vector<double> p, std::vector<double> &a)
{
    double pi=3.141592654;
    /*Input is a vector of PARSEC parameters p=[p1, p2, ...pn] where
    p1=rle
    p2=Xup
    p3=Yup
    p4=YXXup
    p5=Xlow
    p6=Ylow
    p7=YXXlow
    p8=yte
    p9=delta yte (t.e. thickness)
    p10=alpha te
    p11=beta te*/
    Eigen::VectorXd C1(6);
     C1<< 1,1,1,1,1,1;
    Eigen::VectorXd C2(6);
     C2<< pow(p[1],0.5),pow(p[1],1.5),pow(p[1],2.5),pow(p[1],3.5),pow(p[1],4.5),pow(p[1],5.5);
    Eigen::VectorXd C3(6);
     C3<< 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
    Eigen::VectorXd C4(6);
     C4<< (0.5)*pow(p[1],(-0.5)), (1.5)*pow(p[1],(0.5)),(2.5)*pow(p[1],(1.5)),(3.5)*pow(p[1],(2.5)),(4.5)*pow(p[1],(3.5)),(5.5)*pow(p[1],(4.5)) ;
    Eigen::VectorXd C5(6);
     C5<< (-0.25)*pow(p[1],(-1.5)),(0.75)*pow(p[1],(-0.5)),(3.75)*pow(p[1],(0.5)),(8.75) *pow(p[1],(1.5)),(13.25)*pow(p[1],(2.5)),(24.75)*pow(p[1],(3.5));
    Eigen::VectorXd C6(6);
     C6<< 1,0,0,0,0,0;
    Eigen::MatrixXd Cup(6,6);
    Cup.row(0)=C1;Cup.row(1)=C2;Cup.row(2)=C3;Cup.row(3)=C4;Cup.row(4)=C5;Cup.row(5)=C6;

    Eigen::VectorXd C7(6);
     C7<< 1,1,1,1,1,1;
    Eigen::VectorXd C8(6);
     C8<<  pow(p[4],0.5),pow(p[4],1.5),pow(p[4],2.5),pow(p[4],3.5),pow(p[4],4.5),pow(p[4],5.5);
    Eigen::VectorXd C9(6);
     C9<< 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
    Eigen::VectorXd C10(6);
     C10<< (0.5)*pow(p[4],(-0.5)), (1.5)*pow(p[4],(0.5)),(2.5)*pow(p[4],(1.5)),(3.5)*pow(p[4],(2.5)),(4.5)*pow(p[4],(3.5)),(5.5)*pow(p[4],(4.5)) ;
    Eigen::VectorXd C11(6);
     C11<< (-0.25)*pow(p[4],(-1.5)),(0.75)*pow(p[4],(-0.5)),(3.75)*pow(p[4],(0.5)),(8.75) *pow(p[4],(1.5)),(13.25)*pow(p[4],(2.5)),(24.75)*pow(p[4],(3.5));
    Eigen::VectorXd C12(6);
     C12<< 0,0,0,0,0,1;
    Eigen::MatrixXd Clo(6,6);
     Clo.row(0)=C7;Clo.row(1)=C8;Clo.row(2)=C9;Clo.row(3)=C10;Clo.row(4)=C11;Clo.row(5)=C12;
    Eigen::VectorXd bup(6);
     bup<< p[7]+p[8]/2 , p[2], tan((p[9]-p[10]/2)*(pi/180)), 0 , p[3] , (sqrt(2*p[0]));
    Eigen::VectorXd blo(6);
     blo<< p[7]+p[8]/2 , p[5], tan((p[9]-p[10]/2)*(pi/180)), 0 , p[6] , (sqrt(2*p[0]));
    Eigen::VectorXd aup(6);
    Eigen::VectorXd alower(6);
    aup=Cup.inverse()*bup;
    alower=Clo.inverse()*blo;
    for(int i=0;i<12;i++)
    {
        if(i<6)
        {
            a[i]=aup(i);
        }
        else if(i>=6)
        {
            a[i]=alower(i-6);
        }

    }

}

void MainWindow::Panel_Method_Solver(std::vector<double> P, double AOA,double uinf,int NPanel, double &Cl, double &Cdp,double &max_Thickness)
{
       double pi=3.141592654;
       //AOA input in deg
       double dBeta=pi/NPanel;
       double beta=0;

       //std::vector<double>P={0.0155 ,0.296632, 0.060015, -0.4515, 0.296632, -0.06055, 0.453, 0, 0.001260 ,0, 7.36};
       std::vector<double>a(12);
       AOA=AOA*(pi/180);
       double Uinf=uinf*cos(AOA);
       double Vinf=uinf*sin(AOA);

       int k=floor(pi/dBeta);
       std::vector<double>X0(k+1);
       std::vector<double>X0_reversed(k+1);
       std::vector<double>Z_d0(k+1); //lower coordinates
       std::vector<double>Z_u0(k+1); //upper coordinates
       max_Thickness=0;
       Parsec(P,a);

      for (int i=0;i<(k+1);++i)
       {
           X0[i]=(1-cos(beta))/2;
           X0_reversed[i]=(1-cos(beta))/2;
           Shape_Function_PARSIK(a,X0[i],Z_u0[i],Z_d0[i]);

           double Thickness=abs(Z_u0[i]-Z_d0[i]);
           if (Thickness>max_Thickness) max_Thickness=Thickness;

           beta=beta+dBeta;
          // cout<<i<<"   " <<Z_u0[i]<<"   " <<Z_d0[i]<<"   " <<max_Thickness<<endl;

       }
       reverse(Z_d0.begin(),Z_d0.end());
       reverse(X0_reversed.begin(),X0_reversed.end());

       std::vector<double>big_X(2*k+1);//global X coordinates
       std::vector<double>big_Z(2*k+1);//global Z coordinates
       std::vector<double>dx(2*k);
       std::vector<double>dz(2*k);
       std::vector<double>alpha(2*k);
       std::vector<double>cos_theta(2*k);
       std::vector<double>sin_theta(2*k);
       std::vector<double>big_X2(2*k);//Panel length
       Eigen::MatrixXd RHS(2*k,1);//Right Hand Side

       beta=0;
       int count=0; //counter



       for (int i=0;i<(2*k+1);++i)
       {


           if(i<=k)
           {
               big_X[i]=X0_reversed[i];
               big_Z[i]=Z_d0[i];
           }
           else if(i>k)
           {
               big_X[i]=X0[count+1];
               big_Z[i]=Z_u0[count+1];
               count=count+1;
           }
           if(i>0)
           {
               dx[i-1]=big_X[i]-big_X[i-1];
               dz[i-1]=big_Z[i]-big_Z[i-1];
               alpha[i-1]=atan2(dz[i-1],dx[i-1]);
               cos_theta[i-1]=cos(alpha[i-1]);
               sin_theta[i-1]=sin(alpha[i-1]);
               big_X2[i-1]=sqrt(pow(dx[i-1],2)+pow(dz[i-1],2));
               RHS(i-1,0)=-Uinf*sin_theta[i-1] + Vinf*cos_theta[i-1];
           }

           beta=beta+dBeta;
       }



       std::vector<double>x_c(2*k);
       std::vector<double>z_c(2*k);
       for(int i=0;i<2*k;++i)
       {
            x_c[i]=(dx[i]/2)+big_X[i];
            z_c[i]=(dz[i]/2)+big_Z[i];

       }


       Eigen::MatrixXd aij(2*k+1,2*k+1);
       Eigen::MatrixXd bij(2*k,2*k);
       aij.setZero();
       bij.setZero();
       std::vector<double>XP(2*k);
       std::vector<double>ZP(2*k);
       std::vector<double>thetaj(2*k);
       std::vector<double>thetaj1(2*k);
       std::vector<double>R12(2*k);
       std::vector<double>R22(2*k);
       std::vector<double>f(2*k);
       double theta_1 ;
       double theta_2 ;
       for(int i=0;i<2*k;i++)
       {
           for(int j=0;j<2*k;j++)
           {

                  XP[j] = cos_theta[j]*(x_c[i] - big_X[j]) + sin_theta[j]*(z_c[i] - big_Z[j]);
                  ZP[j] = -sin_theta[j]*(x_c[i] - big_X[j]) + cos_theta[j]*(z_c[i] - big_Z[j]);
                  thetaj[j] = atan2(ZP[j],XP[j]);
                  thetaj1[j] = atan2(ZP[j],XP[j] - big_X2[j]);
                  R12[j]=pow(XP[j],2)+pow(ZP[j],2);
                  R22[j]=pow(XP[j] - big_X2[j],2)+pow(ZP[j],2);
                  f[j]=(XP[j]*log(R12[j]) - (XP[j] - big_X2[j])*log(R22[j]) + 2*ZP[j]*(thetaj1[j] - thetaj[j]));
                  if(j!=i)
                  {
                         aij(i,j)=-1/(2*pi)*(thetaj1[j]-thetaj[j]);
                         bij(i,j)=1/(4*pi )*f[j];

                  }
                 else if(j==i)
                 {
                         aij(i,j)=0.5;
                         bij(i,j) = 1/(2*pi)*XP[j]*log(R12[j]);
                 }

                 theta_1 = atan2((z_c[i]-big_Z[2*k]),(x_c[i]-big_X[2*k]));
                 theta_2 = atan2(z_c[i]-big_Z[2*k],x_c[i]-10000000*big_X[2*k]);


           }

           aij(i,2*k)= 1/(2*pi)*(theta_1 - theta_2);
       }


      ///------------A
       Eigen::MatrixXd A= Eigen::MatrixXd::Zero(2*k,2*k);

       for(int i=0;i<2*k;i++)
       {
           A(i,0)=aij(i,0)-aij(i,2*k);

           A(i,(2*k)-1)=aij(i,(2*k)-1)+aij(i,2*k);

       }
       for(int i=0;i<(2*k);i++)
       {
           for(int j=1;j<(2*k)-1;j++)
           {
               A(i,j)=aij(i,j);
           }

       }


     Eigen::VectorXd RHS_new(2*k);
     RHS_new=-bij*RHS;
     Eigen::VectorXd gamma(2*k);
     Eigen::MatrixXd A_inv=A.inverse();
     gamma=A_inv*RHS_new;

     double Ueinf[2*k];
     double Ue[(2*k)-1];
     double Cp[(2*k)-1];

    for(int i=0;i<2*k;i++)
    {
       Ueinf[i]=Uinf*x_c[i]+z_c[i]*Vinf+gamma(i);

    }
    for (int w=1;w<((2*k));++w)
    {
        Ue[w-1]=2*(Ueinf[w-1]-Ueinf[w])/(big_X2[w-1]+big_X2[w]);
        Cp[w-1]=1- pow(Ue[w-1],2) /(pow(Uinf,2)+pow(Vinf,2));
   //        cout<<w-1 << " "<<(Cp[w-1])<<endl;
    }



    double fx[(2*k)-1];
    double fy[(2*k)-1];
    double mj[(2*k)-1];
    for(int j=1;j<(2*k);++j)
    {
        fx[j]=Cp[j]*(big_Z[j+1]-big_Z[j]);
        fy[j]=Cp[j]*(big_X[j+1]-big_X[j]);

        mj[j]=-fx[j]*(big_Z[j+1]+big_Z[j])/2+fy[j]*(((big_X[j+1]+big_X[j])/2)-(1/4));
        //cout<<j<< "  "<<fx[j]<<endl;
    }

    double Fxx=0;
    double Fyy=0;
    double M=0;
    for(int i=0;i<(2*k)-1;i++)
    {
        Fxx+=fx[i];
        Fyy+=fy[i];
        M+=mj[i];

    }
    double Fx=-Fxx;
    double Fy=-Fyy;

    Cl=-sin(AOA)*Fx+cos(AOA)*Fy;
    Cdp=Fx*cos(AOA)+Fy*sin(AOA);//*/


    //cout<< "Cl " <<Cl<<endl;
    //cout<< "Cdp "<<Cdp<<endl;

}

void MainWindow::RandP(std::vector<double> P, std::vector<double> Range, std::vector<double> (&P1))
{
    //if(Range.size()==P.size())
  //  {

        // srand(time(NULL));
        P1.push_back(rand()/( RAND_MAX/(2*Range[0]))+(P[0]-Range[0]));

        P1.push_back(rand()/(RAND_MAX/ (2*Range[1]))+(P[1]-Range[1]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[2]))+(P[2]-Range[2]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[3]))+(P[3]-Range[3]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[4]))+(P[4]-Range[4]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[5]))+(P[5]-Range[5]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[6]))+(P[6]-Range[6]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[7]))+(P[7]-Range[7]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[8]))+(P[8]-Range[8]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[9]))+(P[9]-Range[9]));
         //srand(time(0));
        P1.push_back(rand()/(RAND_MAX/ (2*Range[10]))+(P[10]-Range[10]));

    //}
 /*  else if(Range.size()==2*P.size()) ////??????????????????????????????????/
    {
        P1[0]=(-Range[0]+Range[1])*rand()+Range[0];
        P1[1]=(-Range[2]+Range[3])*rand()+Range[2];
        P1[2]=(-Range[4]+Range[5])*rand()+Range[4];
        P1[3]=(-Range[6]+Range[7])*rand()+Range[6];
        P1[4]=(-Range[8]+Range[9])*rand()+Range[8];
        P1[5]=(-Range[10]+Range[11])*rand()+Range[10];
        P1[6]=(-Range[12]+Range[13])*rand()+Range[12];
        P1[7]=(-Range[14]+Range[15])*rand()+Range[14];
        P1[8]=(-Range[16]+Range[17])*rand()+Range[16];
        P1[9]=(-Range[18]+Range[19])*rand()+Range[18];
        P1[10]=(-Range[20]+Range[21])*rand()+Range[20];
    }//*/


}

void MainWindow::Genetic(int genNo,std::vector<double> P0,std::vector<double>Range,double uinf,double AOA,double Npanels,double &Cloriginal,double & Cl_fittest,std::vector<double>(&fittest_fittest))
{

    double Cdporiginal,max_original;
    Panel_Method_Solver(P0,AOA,uinf,Npanels,Cloriginal,Cdporiginal,max_original);

    int popsize=50;
    double transprob=0.05;
    double crossprob=0.75;
    double mutprob=0.2;
    std::vector<std::vector<double>>newpop;
    newpop.resize(popsize,std::vector<double>(11));
    std::vector<double>Cl(popsize,0.0);
    std::vector<std::vector<double>>p;
    p.resize(popsize,std::vector<double>(11));
    std::vector<std::vector<double>>pop;
    pop.resize(popsize,std::vector<double>(11));
    std::vector<double>P1(11,0.0);
    double Cl_sum;
    std::vector<double>fi(popsize);
    std::vector<double>fittest(popsize);
    std::map<double,int> fi_index;
    int transnumber=3;//=ceil(transprob*popsize);
    int crossnumber=37;//=ceil(crossprob*popsize);
    int mutatnumber=10;//ceil(mutprob*popsize);
    std::vector<double>fittest_new(transnumber);
    std::vector<int>index_new(transnumber);

       time_t seed;

       time(&seed);
       srand(seed*10);
    for(int k=0;k<genNo;++k)
    {

          //Cl.clear();
          p.clear();
          //Cl.resize(11);
          p.resize(popsize,std::vector<double>(11));
          qDebug()<<p.size();

        if(k==0)
        {
            for(int i=0;i<popsize;i++)
            {

                 //srand (time(0) );
                 // RandP(P0,Range,P1);
                 /*for(int i=0;i<11;i++)
                 {
                 qDebug()<<i <<P1[i];}*/
                 ////---
                 P1[0]=rand()/( RAND_MAX/(2*Range[0]))+(P0[0]-Range[0]);
                 P1[1]=rand()/( RAND_MAX/(2*Range[1]))+(P0[1]-Range[1]);
                 P1[2]=rand()/( RAND_MAX/(2*Range[2]))+(P0[2]-Range[2]);
                 P1[3]=rand()/( RAND_MAX/(2*Range[3]))+(P0[3]-Range[3]);
                 P1[4]=rand()/( RAND_MAX/(2*Range[4]))+(P0[4]-Range[4]);
                 P1[5]=rand()/( RAND_MAX/(2*Range[5]))+(P0[5]-Range[5]);
                 P1[6]=rand()/( RAND_MAX/(2*Range[6]))+(P0[6]-Range[6]);
                 P1[7]=rand()/( RAND_MAX/(2*Range[7]))+(P0[7]-Range[7]);
                 P1[8]=rand()/( RAND_MAX/(2*Range[8]))+(P0[8]-Range[8]);
                 P1[9]=rand()/( RAND_MAX/(2*Range[9]))+(P0[9]-Range[9]);
                 P1[10]=rand()/( RAND_MAX/(2*Range[10]))+(P0[10]-Range[10]);
                 /// ----*/
                 double Clnew,Cdpnew,max_new;
                 Panel_Method_Solver(P1,AOA,uinf,Npanels,Clnew,Cdpnew,max_new);

                 if(max_new>0.1)
                    {
                       Clnew=Cloriginal;
                    }
                 else if(max_new<0.01)
                    {
                        Clnew=Cloriginal;
                    }//*/
                // Cl.push_back(Clnew);
                 Cl[i]=Clnew;
                 p[i]=P1;
                 //p.push_back(P1);
             }

         }
        else if(k>0)
        {
            for(int i=0;i<popsize;++i)
            {
                P1=(newpop[i]);
                double Clnew,Cdpnew,max_new;
                Panel_Method_Solver(P1,AOA,uinf,Npanels,Clnew,Cdpnew,max_new);
                //Cl.push_back(Clnew);

                Cl[i]=Clnew;
                // p.push_back(P1);
               //qDebug()<< k<<i <<Clnew;
                p[i]=P1;

            }


        }
        //pop.clear();
        //pop=p;
        for(int j=0;j<popsize;j++)
               {
               for(int i=0;i<11;++i){

               pop[j][i]=p[j][i];
               //qDebug()<<k <<j<<i  <<pop[j][i];;

               }}

        for(int i=0;i<popsize;++i)
             {
                 if(Cl[i]<=Cloriginal)
                    {
                        Cl[i]=Cloriginal;

                    }


             }
        for(int i=0;i<popsize;++i)
        {
                Cl_sum+=Cl[i];


        }



        for(int i=0;i<popsize;++i)
           {
                fi[i]=Cl[i]/Cl_sum;
                fittest[i]=Cl[i]/Cl_sum;

                fi_index[fi[i]]=i;

           }

        sort(fittest.begin(),fittest.end());

        for(int i=0;i<transnumber;++i)
           {
               fittest_new[i]=fittest[i];
               index_new[i]=fi_index[fittest[i]];

           }

       if(k!=(genNo-1))
        {

            for(int i=0;i<transnumber;i++)
            {
                int x=index_new[i];
                newpop[i]=pop[x];
                // qDebug()<<index_new[i];

            }
            //crossOver--------------/

             for(int i=transnumber;i<(crossnumber+transnumber);++i)
              {
                //srand(time(0));
                int max=11;
                int min=0;
                int range = max - min + 1;
                int num = rand() % range + min;     //the point where chromosomes should exchange the binary strings
                int indv1_index=rand()%(pop[0].size()+1);
                int indv2_index=rand()%(pop[0].size()+1);
                std::vector<double>new_child1(11);
                std::vector<double>new_child2(11);
                std::vector<double>Parent1=pop[indv1_index];
                std::vector<double>Parent2=pop[indv2_index];
                for(int i=0;i<11;++i)
                {
                    if(i<num+1)
                    {
                        new_child1[i]=(Parent1[i]);
                        new_child2[i]=(Parent2[i]);
                    }
                    else
                    {
                        new_child1[i]=(Parent2[i]);
                        new_child2[i]=(Parent1[i]);
                    }
                 }

                newpop[i]=new_child1;
                newpop[i]=new_child2;


             }//*/
             //mutation
             for(int i=(crossnumber+transnumber);i<(crossnumber+transnumber+mutatnumber);++i)
             {

                 int indv_index=rand()%(pop[0].size()+1);
                 int index = rand() % 12;
                 std::vector<double>indv(11);
                 std::vector<double>pmut(11);
                 indv=pop[indv_index];
                 RandP(P0,Range,pmut);
                 indv[index]=pmut[index];
                 newpop[i]=indv;
                 //qDebug()<<i+40;

             }
        }

      for(int j=0;j<newpop.size();j++)
       {
       for(int i=0;i<11;++i){

       qDebug()<<k <<j << i << " "<<newpop[j][i];
       }}//*/
       //qDebug()<<pop.size();

    }


    int fittest_index=fi_index[fittest[0]];

    fittest_fittest.resize(11);
    fittest_fittest=pop[fittest_index];
    Cl_fittest=Cl[fittest_index];
    if(Cl_fittest==Cloriginal)
    {
        fittest_fittest=P0;

    }






}





void MainWindow::on_pushButton_clicked()
{

    double P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11;
    P1=ui->P1->text().toDouble();P4=ui->P4->text().toDouble();P6=ui->P6->text().toDouble();P8=ui->P8->text().toDouble();P10=ui->P10->text().toDouble();
    P2=ui->P2->text().toDouble();P5=ui->P5->text().toDouble();P7=ui->P7->text().toDouble();P9=ui->P9->text().toDouble();P11=ui->P11->text().toDouble();
    P3=ui->P3->text().toDouble();
    std::vector<double>P0={P1 ,P2,P3,P4,P5,P6,P7,P8,P9 ,P10,P11};
    double R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11;
    R1=ui->R1->text().toDouble();R4=ui->R4->text().toDouble();R6=ui->R6->text().toDouble();R8=ui->R8->text().toDouble();R10=ui->R10->text().toDouble();
    R2=ui->R2->text().toDouble();R5=ui->R5->text().toDouble();R7=ui->R7->text().toDouble();R9=ui->R9->text().toDouble();R11=ui->R11->text().toDouble();
    R3=ui->R3->text().toDouble();
    std::vector<double>Range={R1 ,R2,R3,R4,R5,R6,R7,R8,R9 ,R10,R11};
    double uinf=ui->uinf->text().toDouble();
    double AOA=ui->AOA->text().toDouble();


    std::vector<double> fittest_fittest(11);
    double Cloriginal;
    double Cl_fittest;
    double Cl_fittest2;
    for(int i=0;i<2;i++){
    Genetic(2, P0,Range,uinf,AOA,50,Cloriginal,Cl_fittest,fittest_fittest);
    P0=fittest_fittest;

    //Genetic(2, fittest_fittest,Range,uinf,AOA,50,Cl_fittest,Cl_fittest2,fittest_fittest);
    }
    QString response=QString("Cl Original=%1").arg(Cloriginal);
    QString response2=QString("Cl Optimized =%1").arg(Cl_fittest);

    ui->cl->setText(response);
    ui->cl_2->setText(response2);



}

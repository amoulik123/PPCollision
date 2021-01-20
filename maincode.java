import java.util.*;
import java.util.Random;
class HiggsFind
{
    // Mass of the Higgs particle has been assumed to be 125 GeV
    Random random= new Random();
    double px,py,pz,px1,py1,pz1;//The momentum of the Photons in the respective planes
    double phi,theta,netenergy1,phix,thetax,netenergy2;//The angles(azimuth and altitude) between the 2 photons
    double mass=125;//in GeV
    double b1,b2,b3,g1,g2,g3;
    double mat1[][]=new double[4][1];//Matrix to store the Energy and momenta of the First Photon
    double mat1x[][]=new double[4][1];//Matrix to store the Energy and momenta of the Second Photon
    double mat2[][]=new double[4][4];// Matrix to store the Boost along X direction
    double mat3[][]=new double[4][4];// Matrix to store the Boost along Y direction
    double mat4[][]=new double[4][4];// Matrix to store the Boost along Z direction
      
    void EnergyOfRestPhoton()//Module 1
     {
        double E=mass/2;//energy of each photon is same as half of the assumed mass of the Higgs Particle(125 GeV)
        double twopi=4*(Math.asin(1.0));
        
        // Energy of Photon
        
        //Calculation of theta
        double costheta=(2*Math.random()) -1.0 ; //generate uniform random in range(-1,1)
        theta=Math.acos(costheta);// in radians
        
        // Calculation of phi
        phi=(twopi*Math.random());//in radians
                
        // Calculation of Momenta of Photon
        px=E*Math.sin(theta)*Math.cos(phi);//Momentum of Photon in the x direction
        py=E*Math.sin(theta)*Math.sin(phi);//Momentum of Photon in the y direction
        pz=E*Math.cos(theta);//Momentum of Photon in the z direction
        netenergy1=(2*Math.sqrt(px*px+py*py+pz*pz)); 
       
        //Energy and Momenta of 1st Photon     
          mat1[0][0]=E;
          mat1[1][0]=px;
          mat1[2][0]=py;
          mat1[3][0]=pz;
         
        //Energy and Momenta of 2nd Photon        
          mat1x[0][0]=E;
          mat1x[1][0]=-px;
          mat1x[2][0]=-py;
          mat1x[3][0]=-pz; 
    }
    
    void EnergyOfRestHiggs()//Module 2
    {
        double pxh,pyh,pzh;//Momenta of the higgs particle in the respective coordinate planes
        
        //Calculation of the momenta of the particle by random generation
        pxh=(100*Math.random());
        pyh=(100*Math.random());
        pzh=(100*Math.random());
        
        double netMomentumh=Math.sqrt(pxh*pxh+pyh*pyh+pzh*pzh);// of the particle
        double energyh=Math.sqrt(mass*mass+(netMomentumh*netMomentumh));//energy of the Higgs Particle
        double vxh=pxh/energyh;//similarly for others. Velocity of the Higgs particle expressed in natural units in the lab frame
        double vyh=pyh/energyh;
        double vzh=pzh/energyh;
        
        //Calculation of Beta and Gamma Values        
        b1=vxh;
        b2=vyh;
        b3=vzh;
        g1=1/(Math.sqrt(1-(b1*b1)));
        g2=1/(Math.sqrt(1-(b2*b2)));
        g3=1/(Math.sqrt(1-(b3*b3)));
               
        // Assignment of values to matrix X
        mat2[0][0]=    g1; mat2[0][1]=-g1*b1; mat2[0][2]=0; mat2[0][3]=0;
        mat2[1][0]=-g1*b1; mat2[1][1]=    g1; mat2[1][2]=0; mat2[1][3]=0;
        mat2[2][0]=     0; mat2[2][1]=     0; mat2[2][2]=1; mat2[2][3]=0;
        mat2[3][0]=     0; mat2[3][1]=     0; mat2[3][2]=0; mat2[3][3]=1;
        
        //Assignment of values to matrix Y
        mat3[0][0]=    g2; mat3[0][1]=0; mat3[0][2]=-g2*b2; mat3[0][3]=0;
        mat3[1][0]=     0; mat3[1][1]=1; mat3[1][2]=     0; mat3[1][3]=0;
        mat3[2][0]=-g2*b2; mat3[2][1]=0; mat3[2][2]=    g2; mat3[2][3]=0;
        mat3[3][0]=     0; mat3[3][1]=0; mat3[3][2]=     0; mat3[3][3]=1;
        
        //Assignment of values to matrix Z
        mat4[0][0]=g3;mat4[0][1]=0;mat4[0][2]=0;mat4[0][3]=-g3*b3;
        mat4[1][0]=0;mat4[1][1]=1;mat4[1][2]=0;mat4[1][3]=0;
        mat4[2][0]=0;mat4[2][1]=0;mat4[2][2]=1;mat4[2][3]=0;
        mat4[3][0]=-g3*b3;mat4[3][1]=0;mat4[3][2]=0;mat4[3][3]=g3;        
    }
    
    public double[][] multiplyMatrices(double[][] firstMatrix, double[][] secondMatrix, int r1, int c1, int c2) {
        double[][] product = new double[r1][c2];
        for(int i = 0; i < r1; i++) {
            for (int j = 0; j < c2; j++) {
                for (int k = 0; k < c1; k++) {
                    product[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                }
            }
        }
        return product;
    }
    
    void EnergyOfHiggsLab()
    {
        Random r=new Random();
        double product1[][] = new double[4][4];
        product1=multiplyMatrices(mat2,mat3,4,4,4);// X*Y
        
        double product2[][] = new double[4][4];
        product2=multiplyMatrices(mat3,mat4,4,4,4);//X*Y*Z
        
        double product3[][]= new double[4][1];//Energy 1       
        product3=multiplyMatrices(mat4,mat1,4,4,1);
        
        double product4[][]=new double[4][1];//Energy 2
        product4=multiplyMatrices(mat4,mat1x,4,4,1);        
                
        double energynew = product3[0][0];
        double pxnew = product3[1][0];
        double pynew = product3[2][0];
        double pznew = product3[3][0];
        
        double energynew1= product4[0][0];
        double pxnew1 = product4[1][0];
        double pynew1 = product4[2][0];
        double pznew1 = product4[3][0];
        
        System.out.println(" Energy 1 is: "+energynew);
        System.out.println(" px 1 is: "+pxnew);
        System.out.println(" py 1 is: "+pynew);
        System.out.println(" pz 1 is: "+pznew);
        System.out.println();
        System.out.println(" Energy 2 is: "+product4[0][0]);
        System.out.println(" px 2 is: "+product4[1][0]);
        System.out.println(" py 2 is: "+product4[2][0]);
        System.out.println(" pz 2 is: "+product4[3][0]);
        
        //Calculation of costheta12
        double costheta12=(pxnew*pxnew1+pynew*pynew1+pznew*pznew1)/((Math.sqrt(pxnew*pxnew+pynew*pynew+pznew*pznew))*Math.sqrt(pxnew1*pxnew1+pynew1*pynew1+pznew1*pznew1));
        System.out.println("The value of cos theta 12 is: "+costheta12);
               
        double sum=0.0;
        double masses[]=new double[10000];
        //Adding computational errors
        for(int i=0;i<50;i++)
        {
            double error1=r.nextGaussian()*0.03+0;//Changing by 0.03 GeV
            double error2=r.nextGaussian()*0.03+0;
            costheta12=Math.cos((Math.acos(costheta12)+(r.nextGaussian()*0.08+0)));//Changing by 0.08 radians
            double e1=energynew+error1,e2=energynew1+error2;
            double mass1=Math.sqrt(2*e1*e2*(1-costheta12));
            sum+=mass1;
            System.out.println(mass1);
            masses[i]=mass1;                    
        }
         double sd=calculateSD(masses);
        System.out.println("The SD is: "+sd+"  The Avg is: "+(sum/10000));
        System.out.println("The assumed mass of the particle was: "+mass);
    }
    
    public double[][] multiply(double a[][], double b[][])
    {
       int rowsInA = a.length;
       int columnsInA = a[0].length; // same as rows in B
       int columnsInB = b[0].length;
       double[][] c = new double[rowsInA][columnsInB];
       for (int i = 0; i < rowsInA; i++) {
           for (int j = 0; j < columnsInB; j++) {
               for (int k = 0; k < columnsInA; k++) {
                   c[i][j] = c[i][j] + a[i][k] * b[k][j];
               }
           }
       }
       return c;
   }

   public double calculateSD(double numArray[])
   {
        double sum = 0.0, standardDeviation = 0.0;
        int length = numArray.length;
        for(int i=0;i<length;i++) {
            sum += numArray[i];
        }
        double mean = sum/length;
        for(int i=0;i<length;i++) {
            standardDeviation += Math.pow(numArray[i] - mean, 2);
        }
        return Math.sqrt(standardDeviation/length);
   }
   
   void main()
   {
       HiggsFind obj=new HiggsFind();
       obj.EnergyOfRestPhoton();
       obj.EnergyOfRestHiggs(); 
       obj.EnergyOfHiggsLab();
   }
}    

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.FileNotFoundException;

/**
 * 
 * CSC340 - Project 1
 * @author Y M. Liu (Simone)
 *
 */
public class Project1 {
	
	private static ArrayList<Matrix> class1 = new ArrayList<Matrix>();
	private static ArrayList<Matrix> class2 = new ArrayList<Matrix>();

    @SuppressWarnings("resource")
	public static void main(String[] args) throws FileNotFoundException {
        Matrix a1 = new Matrix();
        
        //----------------------Testing------------------------------
        /*double [][] data = new double[][]{{4,-2,1}, {8,5,-4}, {-3,1,5}};
        Matrix m3 = new Matrix(data);
        System.out.println("Testing Matrix m3 is a:");
        System.out.println(m3);    
        
        Matrix b = new Matrix(new double[][] {{11},{14},{10}});
        Matrix testGauss = m3.gaussJordan(b);
        System.out.println("Testing GaussJordan");
        System.out.println(testGauss);
        
        ArrayList<Matrix> w1 = new ArrayList<Matrix>();
        ArrayList<Matrix> w2 = new ArrayList<Matrix>();
        
        Matrix w1_vector = new Matrix(new double[][] {{2, 6}, {3, 4}, {3, 8}, {4, 6}});
        Matrix w2_vector = new Matrix(new double[][] {{1, -2}, {3, 0}, {3, -4}, {5, -2}});
        
       
        w1 = w1_vector.toClass_transpose();            
        Matrix cov_test1 = a1.find_covariance(w1);
        System.out.println("Testing Covariance w1");
        System.out.println(cov_test1);
        
        w2 = w2_vector.toClass_transpose();            
        Matrix cov_test2 = a1.find_covariance(w2);
        System.out.println("Testing Covariance w2");
        System.out.println(cov_test2);
        
        Matrix i1 = cov_test1.find_inverse();
        System.out.println("Testing Inverse cov_test1");
        System.out.println(cov_test1.find_inverse());
        
        Matrix i2 = cov_test2.find_inverse();
        System.out.println("Testing Inverse cov_test2");
        System.out.println(cov_test2.find_inverse());
        
        double d1 = cov_test1.find_determinant();
        System.out.println("Testing Determinant cov_test1");
        System.out.println(cov_test1.find_determinant());
        
        double d2 = cov_test2.find_determinant();
        System.out.println("Testing Determinant cov_test2");
        System.out.println(cov_test2.find_determinant());
        
        System.out.println("");
        
        Matrix me1 = a1.find_mean(w1);  
        Matrix me2 = a1.find_mean(w2); 
        
        System.out.println(a1.find_discriminant(me1,me1,i1,d1));
        System.out.println(a1.find_discriminant(me1,me2,i2,d2));
        System.out.println(a1.find_discriminant(me2,me1,i1,d1));
        System.out.println(a1.find_discriminant(me2,me2,i2,d2));*/
        
        
        //------------------------------------------------------------------------------
        
        //Read in the text file
        FileReader inputFile = new FileReader("Resources/p1data.txt");
		Scanner fileReader = new Scanner(inputFile);
		String line = "";
		
		//While the text file is not empty
		while (fileReader.hasNext()) {
			line = fileReader.nextLine();
			Scanner stringReader = new Scanner(line);
			
			//Create 2x1 vectors
			double[][] vector1 = new double[2][1];
			double[][] vector2 = new double[2][1];
			
			//Add vectors to class1 and class2
			vector1[0][0] = stringReader.nextDouble();
			vector1[1][0] = stringReader.nextDouble();
			class1.add(new Matrix(vector1));
			
			vector2[0][0] = stringReader.nextDouble();
			vector2[1][0] = stringReader.nextDouble();
			class2.add(new Matrix(vector2));
		}

		fileReader.close();	
        
		System.out.println("-----------------------------------------------------------------------------------");
        
        //1. Find the mean vectors m1 and m2 for each of the class1 and class2        
        Matrix mean1 = a1.find_mean(class1);
        System.out.println("Mean vector of class1:");
        System.out.println(mean1);
        
        Matrix mean2 = a1.find_mean(class2);
        System.out.println("Mean vector of class2:");
        System.out.println(mean2);     
        
        
		System.out.println("-----------------------------------------------------------------------------------");
        //2. Find the covariance matrices Z1 and Z2 for the classes  
        //Class 1
        Matrix cov1 = a1.find_covariance(class1);
        System.out.println("Covariance of class1:");
        System.out.println(cov1);
        //Class 2
        Matrix cov2 = a1.find_covariance(class2);
        System.out.println("Covariance of class2:");
        System.out.println(cov2);  
		System.out.println("-----------------------------------------------------------------------------------");
        
        //3. Find the determinants of covariance Z1 and Z2
        double deter1 = cov1.find_determinant();
        System.out.println("Determinant of Z1:");
        System.out.println(deter1);
        System.out.println();
        
        double deter2 = cov2.find_determinant();
        System.out.println("Determinant of Z2:");
        System.out.println(deter2);
        System.out.println();
        
		System.out.println("-----------------------------------------------------------------------------------");
        
        //4. Find the determinants of covariance Z1 and Z2
        Matrix inverse1 = cov1.find_inverse();
        System.out.println("Inverse of covariance Z1:");
        System.out.println(inverse1);
        
        System.out.println("Inverse of covariance Z2:");
        Matrix inverse2 = cov2.find_inverse();
        System.out.println(inverse2);
        
		System.out.println("-----------------------------------------------------------------------------------");
        
		//6.Classify means 1 and 2
		System.out.println("Classify classes for the points m1 and m2");
        System.out.println("m1 in class1");
        Matrix m1class1 = a1.find_discriminant(mean1,mean1,inverse1,deter1);
        System.out.println(m1class1);
        
        System.out.println("m1 in class2");
        Matrix m1class2 = a1.find_discriminant(mean1,mean2,inverse2,deter2);
        System.out.println(m1class2);
        
        System.out.println("m2 in class1");
        Matrix m2class1 = a1.find_discriminant(mean2,mean1,inverse1,deter1);
        System.out.println(m2class1);
        
        System.out.println("m2 in class2");
        Matrix m2class2 = a1.find_discriminant(mean2,mean2,inverse2,deter2);
        System.out.println(m2class2);
		
		System.out.println("-----------------------------------------------------------------------------------");
		
		//7, 8.
		ArrayList<Matrix> b_point1 = new ArrayList<Matrix>();
		ArrayList<Matrix> miss_class1 = new ArrayList<Matrix>();
		for (int i = 0; i<class1.size(); i++) {
			Matrix m1c1 = a1.find_discriminant(class1.get(i),mean1,inverse1,deter1);
			Matrix m1c2 = a1.find_discriminant(class1.get(i),mean2,inverse2,deter2);
			if (m1c1.getValue(0, 0)<m1c2.getValue(0, 0)) {				
				System.out.println("g1: "+m1c1.getValue(0, 0));
				System.out.println("g2: "+m1c2.getValue(0, 0));
				miss_class1.add(class1.get(i));	
			}
			a1.boundary_plot(class1.get(i), m1c1, m1c2, b_point1);				
		}
		System.out.println("------Misclassified points Class 1: "+miss_class1.size());
		for (int i=0; i<miss_class1.size(); i++) {
			System.out.println(miss_class1.get(i));			
		}
		System.out.println("Class 1 Boundary points: "+b_point1.size());
		/*for (int i=0; i<b_point1.size(); i++) {
				System.out.println(b_point1.get(i).getValue(1, 0));			
		}*/
		
		
		ArrayList<Matrix> b_point2 = new ArrayList<Matrix>();
		ArrayList<Matrix> miss_class2 = new ArrayList<Matrix>();
		for (int i = 0; i<class2.size(); i++) {
			Matrix m2c1 = a1.find_discriminant(class2.get(i),mean1,inverse1,deter1);
			Matrix m2c2 = a1.find_discriminant(class2.get(i),mean2,inverse2,deter2);
			if (m2c1.getValue(0, 0)>m2c2.getValue(0, 0)) {
				System.out.println("g1: "+m2c1.getValue(0, 0));
				System.out.println("g2: "+m2c2.getValue(0, 0));
				miss_class2.add(class2.get(i));				
			}
			a1.boundary_plot(class2.get(i), m2c1, m2c2, b_point2);				
		}
		System.out.println("------Misclassified points Class 2: "+miss_class2.size());
		for (int i=0; i<miss_class2.size(); i++) {
			System.out.println(miss_class2.get(i));			
		}
		System.out.println("Class 2 Boundary points: "+b_point2.size());
		/*for (int i=0; i<b_point2.size(); i++) {
			System.out.println(b_point2.get(i).getValue(1, 0));			
		}*/
		
		System.out.println("-----------------------------------------------------------------------------------");

        //9.
		System.out.println("Estimate solution for linear system:");
        Matrix linear = new Matrix(new double [][]
            	{{2,1,-1,-1,1,0,-1,-1},
            	{1,0,2,0,-1,-2,2,2},
            	{0,-2,5,4,-1,0,3,1},
            	{1,1,-7,3,2,1,-1,0},
            	{1,1,2,3,-2,2,2,9},
            	{0,-3,-2,2,0,2,4,-5,-3},
            	{-2,5,-1,1,1,3,0,-2,4},
            	{1,0,1,1,0,2,1,1,-4}});		

        
        Matrix linear_b = new Matrix(new double [][] {{1},{-1},{2},{-2},{3},{-3},{4},{-4}});    
        System.out.println(linear);
        
        System.out.println("Solution, order [x, y, z, w, a, b, c, d]");
        Matrix solve_linear = linear.gaussJordan(linear_b);
        System.out.println(solve_linear.getCol(8));      
        
        
        double linear_deter = linear.find_determinant();
        System.out.println("Determinant of matrix:");
        System.out.println(linear_deter);
        
        Matrix linear_inverse = linear.find_inverse();
        System.out.println("Inverse of matrix:");
        System.out.println(linear_inverse);
        
        double inverse_deter = linear_inverse.find_determinant();
        System.out.println("Determinant of inverse matrix:");
        System.out.println(inverse_deter);
        
        double deter_product = linear_deter * inverse_deter;
        System.out.println("Product of determinants matrix and inverse matrix");
        System.out.println(deter_product);
        System.out.println("");

        
        Matrix checking = linear.mult(linear_inverse);
        System.out.println("Product of matrix and its inverse");
        System.out.println(checking);
        
		System.out.println("-----------------------------------------------------------------------------------");

		//10.
		System.out.println(linear.find_condition(linear_inverse));
        
     
    }  
    
    //----------------------Functions----------------------------
    

}
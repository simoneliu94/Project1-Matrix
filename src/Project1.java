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

        double [][] data = new double[][]{{4,-2,1}, {8,5,-4}, {-3,1,5}};
        Matrix m3 = new Matrix(data);
        System.out.println("Matrix m3 is a:");
        System.out.println(m3);    
        
        Matrix b = new Matrix(new double[][] {{11},{14},{10}});
        Matrix testGauss = m3.gaussJordan(b);
        System.out.println("Testing GaussJordan");
        System.out.println(testGauss);
        
        
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
        
        //9.
		System.out.println("Estimate solution for linear system:");
        double [][]linear_data = new double [][]
            	{{2,1,-1,-1,1,0,-1,-1},
            	{1,0,2,0,-1,-2,2,2},
            	{0,-2,5,4,-1,0,3,1},
            	{1,1,-7,3,2,1,-1,0},
            	{1,1,2,3,-2,2,2,9},
            	{0,-3,-2,2,0,2,4,-5,-3},
            	{-2,5,-1,1,1,3,0,-2,4},
            	{1,0,1,1,0,2,1,1,-4}};
            Matrix linear = new Matrix(linear_data);
            System.out.println(linear);
            Matrix linear_b = new Matrix(new double [][] {{1},{-1},{2},{-2},{3},{-3},{4},{-4}});
            
            Matrix solve_linear = linear.gaussJordan(linear_b);
            System.out.println(solve_linear);
     
    }  
}
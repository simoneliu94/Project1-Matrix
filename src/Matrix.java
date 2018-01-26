import java.util.ArrayList;
import java.util.Arrays; 

/**
 * 
 * A class of all matrix operations
 * @author Y M. Liu (Simone)
 *
 */
public class Matrix {
    private int numRows = 0;
    private int numCols = 0;
    private double matrix[][];
    
/**
 * 
 */
    public Matrix() {
    	
    }
/**
 * Represents a matrix in 2 dimensions
 * @param matrix this is a 2D array that takes double
 */
    public Matrix(double[][] data)
    {
        this.numRows = data.length;
        this.numCols = data[0].length;
        this.matrix = new double[numRows][numCols];
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                this.matrix[i][j] = data[i][j];
            }
        }        
    }
/**
 * Gets the number of rows and columns and create a new matrix 
 * @param rows is the number of row in the matrix
 * @param cols is the number of column in the matrix
 */
    public Matrix(int rows, int cols)
    {
        this.numRows = rows;
        this.numCols = cols;
        this.matrix = new double[numRows][numCols];
    }
    
    
    public double getValue(int row, int col) {
    	return this.matrix[row][col];
    }
    
    public Matrix getRow(int row) {
    	Matrix A = new Matrix(1, numCols);
    	for (int i=0; i<A.numCols; i++)
        {            
            A.matrix[0][i] = this.matrix[row][i];
        }
    	return A;
    }
    
    public Matrix getCol(int col) {
    	Matrix A = new Matrix(numRows, 1);
    	for (int i=0; i<A.numRows; i++)
        {            
            A.matrix[i][0] = this.matrix[i][col];
        }
    	return A;
    }

    
    public ArrayList<Matrix> toClass(){
    	ArrayList<Matrix> w1 = new ArrayList<Matrix>();
    	for (int i = 0; i<this.numRows; i++) { 
    		w1.add(this.getRow(i));    		
    	}
    	return w1;
    }
    
    public ArrayList<Matrix> toClass_transpose(){
    	Matrix trans = this.transpose();
    	ArrayList<Matrix> w1 = new ArrayList<Matrix>();
    	for (int i = 0; i<trans.numCols; i++) { 
    		w1.add(trans.getCol(i));    		
    	}
    	return w1;
    }

/**
 * Represents a matrix in 2D in String format
 */
    public String toString()
    {
        String output = "";
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                output += matrix[i][j] + " ";                
            }
            output += "\n";
        }

        System.out.println(numRows +"x"+ numCols + " matrix");
        return output;
    }
    
    
/**
 * Swaps rows and columns to each other
 * @return a new matrix after swapping
 */
    public Matrix transpose()
    {
        Matrix A = new Matrix(numCols,numRows);
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                A.matrix[j][i] = this.matrix[i][j];
            }
        }
        return A;
    }    
/**
 * Adds 2 matrices 
 * @param B is a matrix you want to add to
 * @return a matrix result after adding
 */
    public Matrix add(Matrix B)
    {
        Matrix A = this;
        if(B.numRows != A.numRows || B.numCols != A.numCols)
        {
            System.out.println("CANNOT ADD");
            return null;
        }

        Matrix C = new Matrix(numRows,numCols);
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];                
            }
        }
        return C;
    }
 
/**
 * 
 * @param B
 * @return
 */
    public Matrix subtract(Matrix B) {
    	Matrix A = this;
        if(B.numRows != A.numRows || B.numCols != A.numCols)
        {
            System.out.println("CANNOT SUBTRACT");
            return null;
        }

        Matrix C = new Matrix(numRows,numCols);
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];                
            }
        }
        return C;
    }
/**
 * Multiplies 2 matrices    
 * @param B is a matrix you want to multiply to
 * @return a matrix result after multiplying
 */
    public Matrix mult(Matrix B)
    {
        Matrix A = this;
        if (A.numCols != B.numRows)
        {
            System.out.println("CANNOT MULTIPLY");
            return null;
        }
        
        Matrix C = new Matrix(A.numRows,B.numCols);
        for (int i=0; i< C.numRows; i++)
        {
            for (int j=0; j<C.numCols; j++)
            {
                for (int m=0; m<A.numCols; m++)
                C.matrix[i][j] += (A.matrix[i][m] * B.matrix[m][j]);
            }
        }
        return C;        
    }
/**
 * Multiplies every number in the matrix by a same number
 * @param scalar is a number
 * @return a matrix result after multiplying
 */
    public Matrix mult(double scalar)
    {
        Matrix A = this;
        Matrix C = new Matrix(numRows,numCols);
        
        for (int i=0; i<numRows; i++)
        {
            for (int j=0; j<numCols; j++)
            {
                C.matrix[i][j] = A.matrix[i][j]*scalar;             
            }
        }
        return C;
    }
/**
 * Compares 2 matrices    
 * @param B is a matrix you want to compare to
 * @return a result if the matrices are equal or not
 */
    public boolean equals(Matrix B)
    {
        Matrix A = this;
        double firstOne[][] = new double[A.numRows][A.numCols];
        double secondOne[][] = new double[B.numRows][B.numCols];
        return Arrays.deepEquals(firstOne, secondOne);
    }
    
/**
 * Create a square identity matrix    
 * @param size
 * @return
 */
    public Matrix create_identity (int size) {
    	Matrix A = new Matrix(size,size);
    	for (int i = 0; i < size; i++) {
    		for (int j = 0; j < size; j++) {
    			if(i==j) {
    				A.matrix[i][j]=1.0;
    			}
    			else {
    				A.matrix[i][j]=0.0;
    			}
    		}
    	}
    	return A;
    }
    
/**
 * Augment coefficient matrix to create a nx2n or nxn+1 matrix     
 * @param A: given
 * @param I: identity matrix
 * @return
 */
    public Matrix augment(Matrix A, Matrix I) {
    	Matrix C = new Matrix(A.numRows, A.numCols+I.numCols);
    	for (int i = 0; i < C.numRows; i++) {
    		for (int j = 0; j < C.numCols; j++) {
    			if(j < A.numCols) {
    				C.matrix[i][j] = A.matrix[i][j];
    			}
    			else {
    				C.matrix[i][j] = I.matrix[i][j-A.numCols];
    			}
    		}
    	}    	
		return C;    	
    }
 
/**
 * Find the mean vector of a class
 * @param aClass/matrix
 * @return Matrix (the mean of a class) 
 */
    public Matrix find_mean (ArrayList<Matrix> aClass){
    	//Get the first vector out of the class
    	Matrix m1 = aClass.get(0);
    	//Find sum of the class
    	for (int i = 1; i < aClass.size(); i++) {
        	m1 = m1.add(aClass.get(i));
        }       
    	//Find average of the class
        m1 = m1.mult(1.0/aClass.size());
    	
        //Return the mean vector
		return m1;	    	
    }
    
/**
 * Find the covariance matrix
 * @param aClass
 * @return
 */
    public Matrix find_covariance(ArrayList<Matrix> aClass) {
    	//Create a new class to store the result after found the nxn product
    	ArrayList<Matrix> class1 = new ArrayList<Matrix>();
    	
    	//Find mean
    	Matrix mean = find_mean(aClass);
    	
    	//1.Subtract the mean from each vector of the class
    	//2.Multiply the result from 1. to its transpose to find nxn product
    	for (int i = 0; i < aClass.size(); i++) {
    		Matrix m1 = aClass.get(i).subtract(mean);    		
    		Matrix square = m1.mult(m1.transpose());
    		
    		class1.add(square);    		
    	}
    	    	
    	Matrix covariance = find_mean(class1);
    	
		return covariance;    	
    }

/**
 * Swaps 2 rows in matrix
 * @param rowA
 * @param rowB
 */
    public void interchange_row(int rowA, int rowB) {
    	double[] temp_row = matrix[rowA];
    	matrix[rowA] = matrix[rowB];
		matrix[rowB] = temp_row;
    }

/**
 * Implement Gauss Jordan elimination    
 * @param b
 * @return
 */
    public Matrix gaussJordan(Matrix b) {
    	@SuppressWarnings("unused")
		int E = 1;   
    	int p;
    	Matrix C = augment(this,b);
    	    	
    	for(int j = 0; j < C.numRows; j++) {
    		p = j;
    		//Compute pivot, find max of the column
    		for(int i = j; i < C.numRows; i++) {
    			if (Math.abs(C.matrix[p][j]) < Math.abs(C.matrix[i][j])) {
    				p = i;
    			}
    		}    	
    		
    		if (C.matrix[p][j] == 0) {
				E = 0;
			}
    		
    		//Swap 2 rows j and p
    		if(p > j) {
    			C.interchange_row(j, p);
    		}
    		
    		//Divide row j by Cjj
    		double Cjj = C.matrix[j][j];    
    		for (int i = 0; i < C.numCols; i++) {    						
				C.matrix[j][i] = C.matrix[j][i] / Cjj;
    		}
    			
    		//For each i != j, subtract Cij times row j from row i
    		for (int i = 0; i < C.numRows; i++) {
    			if (i!=j) {
    				double Cij = C.matrix[i][j];
    				for (int m=0; m < C.numCols; m++) {
    					C.matrix[i][m] = C.matrix[i][m] - (C.matrix[j][m] * Cij);
    				}
    			}					
			} 	 
    	}    	
		return C;    	 
    }

/**
 * Find determinant matrix
 * @return
 */
    public double find_determinant() {
    	int r = 0;
    	int p;
    	double determinant = 1.0;
    	Matrix A = new Matrix(this.matrix);
    	
    	//Double check on finding determinant matrix
    	//double deter = A.matrix[0][0]*A.matrix[1][1] - A.matrix[0][1]*A.matrix[1][0];
    	//System.out.println(deter);
    	    	
    	for(int j = 0; j < A.numRows; j++) {
    		p = j;
    		//Compute pivot, find max of the column
    		for(int i = j; i < A.numRows; i++) {
    			if (Math.abs(A.matrix[p][j]) < Math.abs(A.matrix[i][j])) {
    				p = i;
    			}
    		}   		
    		
    		if (A.matrix[p][j] == 0) {
    			determinant = 0.0;
			}
    		//Swap 2 rows j and p
    		if(p > j) {
    			A.interchange_row(j, p);
    			r++;    			
    		}
    		
    		//subtract Aij /Ajj times row j from row i
    		for (int i = 0; i < A.numRows; i++) {
    			if(i>j) {
    				double Aij = A.matrix[i][j];
    				double Ajj = A.matrix[j][j];
    				for (int m = 1; m < A.numCols; m++) {
    					A.matrix[i][m] = A.matrix[i][m] - (A.matrix[j][m] * (Aij/Ajj));
    				}    				
    			}
    		}
    	}    	
    	//multiply diagonally 
    	for (int i = 0; i < A.numRows; i++) {
			determinant = determinant * A.matrix[i][i];
		}    	
    	
    	//(-1)^r( A11xA22x…xAnn)
    	determinant = determinant * Math.pow((-1), r);

		return determinant;        	
    }
    

/**
 * Find inverse of a matrix    
 * @return
 */
    public Matrix find_inverse() {
		int E = 1;   
    	int p;
    	Matrix C = augment(this,create_identity(numRows));
    	    	
    	for(int j = 0; j < C.numRows; j++) {
    		p = j;
    		//Compute pivot, find max of the column
    		for(int i = j; i < C.numRows; i++) {
    			if (Math.abs(C.matrix[p][j]) < Math.abs(C.matrix[i][j])) {
    				p = i;
    			}
    		}    	
    		
    		if (C.matrix[p][j] == 0) {
				E = 0;
			}
    		
    		//Swap 2 rows j and p
    		if(p > j) {
    			C.interchange_row(j, p);
    		}
    		
    		//Divide row j by Cjj
    		double Cjj = C.matrix[j][j];    
    		for (int i = 0; i < C.numCols; i++) {    						
				C.matrix[j][i] = C.matrix[j][i] / Cjj;
    		}
    			
    		//For each i != j, subtract Cij times row j from row i
    		for (int i = 0; i < C.numRows; i++) {
    			if (i!=j) {
    				double Cij = C.matrix[i][j];
    				for (int m=0; m < C.numCols; m++) {
    					C.matrix[i][m] = C.matrix[i][m] - (C.matrix[j][m] * Cij);
    				}
    			}					
			} 	 
    	} 
    	
    	//Save inverse matrix since the form is now [A,I]
    	Matrix inverse = new Matrix(C.numRows,C.numCols/2);
    	for (int i = 0; i < inverse.numRows; i++) {
    		for (int j = C.numCols/2; j < C.numCols; j++) {
    			inverse.matrix[i][j-C.numCols/2] = C.matrix[i][j];
    		}
    	}    	
		return inverse;    	
    }

/**
 *     
 * @param points
 * @param mean
 * @param inverse
 * @param det
 * @return
 */
    public Matrix find_discriminant(Matrix points, Matrix mean, Matrix inverse, double det) {
    	Matrix g = new Matrix(); 

    	g = points.transpose().subtract(mean.transpose());
    	g = g.mult(-0.5);
    	g = g.mult(inverse);
    	g = g.mult(points.subtract(mean));
    	
    	double g2 = (-0.5*Math.log(det));
    	g2 = g2 + Math.log(0.5);
    	
    	Matrix g3 = new Matrix(new double[][] {{g2}});
    	
    	g = g.add(g3);
    	return g;
    }

/**
 *     
 * @param aClass
 * @param dis1
 * @param dis2
 * @return
 */
    public void boundary_plot(Matrix point, Matrix dis1, Matrix dis2, ArrayList<Matrix> b_point) {
    	//ArrayList<Matrix> plot_points = new ArrayList<Matrix>();
    	double step = 0.05;
    	double eps = 0.1;
    	
    	double g1 = dis1.getValue(0,0)*step;
    	double g2 = dis2.getValue(0,0)*step;
    	double mag = Math.abs(g1-g2);
    	
    	if(mag<eps) {
    		b_point.add(point);
    	}    	
    }    
    
}

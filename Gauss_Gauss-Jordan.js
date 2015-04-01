/*	Global Variables:

x= x Matrix 
numSources= Number of sources
numDest= Number of destinations
*/

// CHECKOUT FUNCTIONALITY

/*
	function: calculate_U_V
	Description: 
		This is the main function of this file. It calculates the dual vectors U and V as result of Gauss-Jordan algorithm. This function is divided in 4 steps:
		1) Data initialization;
		2) Creation and initialization of a square matrix defined in the dual problem;
		3) Gauss algorithm: at first the matrix will be extended by a constant vector, containing the distances between sources and destinations.
		   Than Gauss algorithm will be applied, where the whole system is included in the extended matrix.
		4) Gauss-Jordan elimination through pivoting.
	Parameters:
	Return:
		- U_V [array (numDest + numSources-1)]: Array of  numDest + numSources elements, that contains both vectors U and V as result of Gauss-Jordan algorithm.
*/
function calculate_U_V(){
	DEBUG_MODE=false;
	
//step 1: data initialization.
	var constant_column=new Array(numSources+numDest-1);
	var index=0;	
	var U_V=new Array(numDest + numSources-1);
	var U_V_temp=new Array(numDest + numSources-1);
	
	//Matrix creation and initialization to zero
	var gauss_Matrix=new Array(numSources+numDest-1);
	for (var i=0; i<gauss_Matrix.length; i++)	gauss_Matrix[i]=new Array(numSources+numDest-1);
	
	for (var i=0;i<gauss_Matrix.length; i++){
		for (var j=0; j<gauss_Matrix[i].length; j++)   {
			gauss_Matrix[i][j]=0;
		}	
	}
	
//step 2: Changing gauss_Matrix values in order to apply Gauss algorithm. In 'change_Values' function 'constant_column' values are properly setted.
	for (var i=0;i<numSources;i++){
		for (var j=0;j<numDest;j++){
			if (x[i][j]!=null){
				change_Values(gauss_Matrix,index,i,j,constant_column);
				index+=1;
			}
		}
	}
	if (DEBUG_MODE==true)	alert(gauss_Matrix);
	
//step 3: Gauss algorithm.
	//matrix extension. 
	var gauss_Matrix_Extended=extend_matrix(gauss_Matrix,constant_column);
	if (DEBUG_MODE==true)	alert(gauss_Matrix_Extended);
	// Gauss algorithm
	gauss_Matrix_Extended=gauss_Algorithm(gauss_Matrix_Extended);
	if (DEBUG_MODE==true)	alert(gauss_Matrix_Extended);

//step 4: Gauss-Jordan elimination algorithm.
	U_V=gauss_Jordan_Algorithm(gauss_Matrix_Extended,U_V);
	
	return U_V;
}

/*
	function: gauss_Algorithm
	Description: 
		Implementation of Gauss algorithm. The algorithm is included in a while loop. Logical steps:
		1) If the values of the first column are all zeros, go to step 3. If the first element of the first row is zero, change the row with another one,
		   where the first element is not zero;
		2) If the first element of the first row is not zero, perform row operations contained in 'calculate_Gauss' function;
		3) Now the first column has all zero values, except perhaps the first element. Go to step 1 and consider the submatrix obtained by eliminating the first row and column.
	Parameters:
		- gauss_Matrix_Extended [matrix ((numDest + numSources-1)*(numDest + numSources))]:	system before Gauss algorithm.
	Return:
		- final_Gauss_Matrix [matrix (numDest + numSources-1)*(numDest + numSources))]: 	system after Gauss algorithm. 
*/

function gauss_Algorithm(gauss_Matrix_Extended) {
//data initialization
	var final_Gauss_Matrix=gauss_Matrix_Extended.slice();
	var end_Algorithm=false;
	
//algorithm
	if (gauss_Matrix_Extended.length==1)
		end_Algorithm=true;
	while (end_Algorithm==false){
		var column_to_check=get_First_Column(gauss_Matrix_Extended);
		all_zero_column=check_Column(column_to_check);
		// If the first element of all rows are zero, copy the first row and column of 'gauss_Matrix_Extended' in 'final_Gauss_Matrix'.
		if (all_zero_column==true) {
			final_Gauss_Matrix=write_final_matrix(gauss_Matrix_Extended,final_Gauss_Matrix);
			gauss_Matrix_Extended=calculate_Submatrix(gauss_Matrix_Extended);
		}				
			else{
				if (gauss_Matrix_Extended[0][0]==0){
					gauss_Matrix_Extended=change_Row(gauss_Matrix_Extended,ordered_variables);		
				}
				var column=new Array(gauss_Matrix_Extended.length);
				column=get_First_Column(gauss_Matrix_Extended);
				
				// if the first column does not have zero values except the first one, perform row operations. 
				// Otherwise copy the first row and column in 'final_Gauss_Matrix'.
				var column_NotZero=column_Test_NotZero(column);
				if (column_NotZero==true){
					gauss_Matrix_Extended=calculate_Gauss(gauss_Matrix_Extended);				
				}
				else {
					final_Gauss_Matrix=write_final_matrix(gauss_Matrix_Extended,final_Gauss_Matrix);
					gauss_Matrix_Extended=calculate_Submatrix(gauss_Matrix_Extended);
				}
			}
		if (gauss_Matrix_Extended.length==1){
			final_Gauss_Matrix=write_final_matrix(gauss_Matrix_Extended,final_Gauss_Matrix);
			end_Algorithm=true;
		}
	}
	return final_Gauss_Matrix;
}


/*
	function: gauss_Jordan_Algorithm
	Description: 
		Gauss-Jordan elimination through pivoting
	Parameters:
		- gauss_Matrix_Extended [matrix ((numDest + numSources-1)*(numDest + numSources))]
		- U_V [array (numDest + numSources-1)]								
	Return:
		- U_V [array (numDest + numSources-1)]:	dual variables as result of Gauss-Jordan algorithm.
*/
function gauss_Jordan_Algorithm(gauss_Matrix_Extended,U_V){
	
	var n=gauss_Matrix_Extended.length;
	U_V[U_V.length-1]=gauss_Matrix_Extended[n-1][n]/gauss_Matrix_Extended[n-1][n-1];
	for (var i=n-2; i>=0; i--){
		for (var j=n-1; j>i; j--){
			if (gauss_Matrix_Extended[i][j]!=0){
				gauss_Matrix_Extended[i][n]=gauss_Matrix_Extended[i][n]-(gauss_Matrix_Extended[i][j]*U_V[j]);
			}
		}
		U_V[i]=gauss_Matrix_Extended[i][n]/gauss_Matrix_Extended[i][i];
	}
	return U_V;
}

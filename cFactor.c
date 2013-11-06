#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h> 
 
#define TRUE 1
#define FALSE 0
// To be used in pseudo random generator function f_rho
#define RHO 1 

/* Will hold the factors of number to be factored,
   reset between each run. */
mpz_t* factors;
int nrOfFactors;
int primesCleared;
 
/*
    Adds a factor to factors.

    Args:
        @factor - Factor to be added to factors array.

    Returns: -

*/
void addToFactors(mpz_t factor){
 
    mpz_set(factors[nrOfFactors], factor);
    nrOfFactors++;
 
}
 
/*
    Performs trial division on compositeNumber.
    Will divide out factors found.

    Args: 
        @compositeNumber - Integer to be factored.

    Returns: - 

*/
void trialDivision(mpz_t compositeNumber){
 
    mpz_t test_prime;
    mpz_t upper_bound;
    mpz_t mod_test_prime;
 
    mpz_init(test_prime);
    mpz_init(upper_bound);
    mpz_init(mod_test_prime);
 
    mpz_set_ui(test_prime, 2);
    mpz_set_ui(upper_bound, 1000);
    
    // Stop when current prime to be tested becomes greater than upper_bound 
    while(mpz_cmp(test_prime, upper_bound) <= 0){
        
        // If compositeNumber is a prime, there is no need to continue factorizing
        if(mpz_probab_prime_p(compositeNumber, 2)){
            return;
        }
 
        // Check if test_prime is a factor of compositeNumber
        mpz_mod(mod_test_prime, compositeNumber, test_prime);
        // If it was a factor
        if(mpz_cmp_ui(mod_test_prime, 0) == 0){
            addToFactors(test_prime);
            // Remove the found factor from compositeNumber
            mpz_fdiv_q(compositeNumber, compositeNumber, test_prime);
        }
        else{
            // If it wasn't a factor, time to move on to the next prime
            mpz_nextprime(test_prime, test_prime);
        }
    }
 
}
 
/*

    Pseudo random function used in Pollard rho.
    Performs the function f(x) = x^2 + a mod n.
    Overwrites x with the result.

    Args: 
        @x - Used as in function above, will be overwritten by result.
        @n - Used as in function above.
        @a - Used as in function above.

    Returns: -

*/
void f_rho(mpz_t x, mpz_t n, unsigned int a){
 
    mpz_pow_ui(x, x, 2);
    mpz_add_ui(x, x, a);
    mpz_mod(x, x, n);
 
}
 
/*

    Exponential factoring algorithm. Tries to find a factor
    of compositeNumber. If it successfully finds a factor,
    the factor will be overwritten into res. Otherwise
    res will be overwritten by 0.
    Uses Floyd's cycle detection.

    Args:
        @compositeNumber - Integer to be factored.
        @res - Holder for result.

    Returns: -

*/
void pollardRhoStd(mpz_t compositeNumber, mpz_t res){
 
    unsigned int j = 0;

    // The tortoise
    mpz_t x;
    // The hare
    mpz_t y;
    mpz_t d;
    
    mpz_init(x);
    mpz_init(y);
    mpz_init(d);
 
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
 
    while(TRUE){
        // Advance tortoise one step, hare two steps
        f_rho(x, compositeNumber, 1);
        f_rho(y, compositeNumber, 1);
        f_rho(y, compositeNumber, 1);
    
        // Check if a cycle has been found
        mpz_sub(d, x, y);
        mpz_abs(d, d);
        mpz_gcd(d, d, compositeNumber);
 
        // Cycle was found (-> factor was probably found)
        if(mpz_cmp_ui(d, 1) != 0){
            // Make sure that the factor doesn't equal compositeNumber
            if(mpz_cmp(d, compositeNumber) != 0){
                mpz_set(res, d);
                return;
            }
            // If the factor was equal to compositeNumber,
            // the algorithm has failed.
            else{
                mpz_set_ui(res, 0);
                return;
            }
        }
        
    }
    
    // A cycle was not found in reasonable time
    mpz_set_ui(res, 0);
    return;
 
}

/*
    Puts the min(a,b) into res.
    If a equals b, res will be set to b.

    Args:
        @a - to be compared with b
        @b - to be compared with a
        @res - to be overwritten bi min(a,b)

    Returns: -

*/
void min(mpz_t a, mpz_t b, mpz_t res){

    if(mpz_cmp(a, b) < 0)
        mpz_set(res, a);
    else
        mpz_set(res, b);

}

/*

    Exponential factoring algorithm. Tries to find a factor
    of compositeNumber. If it successfully finds a factor,
    the factor will be overwritten into res. Otherwise
    res will be overwritten by 0.
    Uses Brent's cycle detection.

    Args:
        @compositeNumber - Integer to be factored.
        @res - Holder for result.

    Returns: -
    
*/
void pollardRhoBrent(mpz_t compositeNumber, mpz_t res){

    // Tortoise
    mpz_t x;
    // Hare
    mpz_t y;
    mpz_t d;
    mpz_t r;
    mpz_t m;
    mpz_t q;
    mpz_t ys;
    mpz_t minVal;    
    mpz_t k;
    mpz_t temp;

    mpz_init(x);
    mpz_init(y);
    mpz_init(d);
    mpz_init(r);
    mpz_init(m);
    mpz_init(k);
    mpz_init(q);
    mpz_init(ys);
    mpz_init(minVal);
    mpz_init(temp);

    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
    mpz_set_ui(r, 1);
    mpz_set_ui(m, 100);
    mpz_set_ui(q, 1);
    mpz_set_ui(ys, 1);
    mpz_set_ui(temp, 1);
    unsigned int i;

    while(mpz_cmp_ui(d, 1) == 0){
        // Teleport the tortoise!
        mpz_set(x, y);
        
        // Advance the hare r steps
        for(i = 0; mpz_cmp_ui(r, i) > 0; i++){
            f_rho(y, compositeNumber, RHO);     
        }
        
        mpz_set_ui(k, 0);
       
        // While r is less than k AND the gcd did not yeild a factor
        while((mpz_cmp(k, r) < 0) && (mpz_cmp_ui(d, 1) == 0)){
            mpz_set(ys, y);
            mpz_sub(minVal, r, k);
            min(m, minVal, minVal);
            
            // Try to decrease the number of gcd calls
            for(i = 0; mpz_cmp_ui(minVal, i) > 0; i++){
                f_rho(y, compositeNumber, RHO);
                mpz_sub(temp, x, y);
                mpz_abs(temp, temp);
                mpz_mul(q, temp, q);
                mpz_mod(q, q, compositeNumber);
            }
            
            mpz_gcd(d, q, compositeNumber);
            mpz_add(k, k, m);
        }
        
        mpz_mul_ui(r, r, 2);
    }

    // The factor was equally to compositeNumber, meaning pollard
    // brent failed, time to do regular floyd pollard
    if(mpz_cmp(r, compositeNumber) == 0){
        while(TRUE){
            f_rho(ys, compositeNumber, RHO);
            mpz_sub(temp, x, ys);
            mpz_abs(temp, temp);
            mpz_gcd(d, temp, compositeNumber);
            if(mpz_cmp_ui(d, 1) != 0)
                break;    
        }
    }
      
    // Floyd pollard also failed
    if(mpz_cmp(r, compositeNumber) == 0)
        mpz_set(res, 0);
    // A factor was found
    else
        mpz_set(res, d);
    
    return;
}

/*
    This function tries to find all the prime factors of compositeNumber.

    Args:
        @compositeNumber - The integer to be factored.
    
    Returns:
        TRUE - If all the factors were (probably) found.
        FALSE - If it failed to find all of the factors.

*/
int factorThis(mpz_t compositeNumber){
    
    trialDivision(compositeNumber);
    mpz_t res;
    mpz_init(res);
    while(TRUE){
        if(mpz_probab_prime_p(compositeNumber, 2)){
            addToFactors(compositeNumber);
            return TRUE;
        }
        
        // Choose std for Floyd, Brent for Brent
        pollardRhoStd(compositeNumber, res);
        //pollardRhoBrent(compositeNumber, res);     

        // If res is 0 the Pollard algorithm failed
        if(mpz_cmp_ui(res, 0) == 0){
            return FALSE;    
        }
 
        // Check if the factor provided by Pollard actually is a prime
        if(mpz_probab_prime_p(res, 2)){
            addToFactors(res);    
        }
        else
            return FALSE;
 
        // Divide out found factor
        mpz_fdiv_q(compositeNumber, compositeNumber, res);
    }
    return FALSE;
}
 
/*
    
    Initialize an array to be filled with mpz_t.

    Args:
        array - Array to be initialized.
        size - Size of array to be initialized.

    Returns: -

*/
void initArray(mpz_t* array, int size){
 
    int i;
    for(i = 0; i < size; i++){
        mpz_init(array[i]);
    }
 
}
 
/*

    Print factors found in global factors array.
    Will print as many as indiciated by nrOfFactors.

    Args: -
    Returns: -

*/
void printFactors(){
    
    int i;
    for(i = 0; i < nrOfFactors; i++){
        gmp_printf("%Zd\n", factors[i]);
    }
    printf("\n");
}

/*
    
    This function is called when main receives a SIGINT.
    Prints out fail on remaining numbers queued to be factored.
    Terminates program.
    
    Args: 
        sig: Signal received.

    Returns: -

*/
void interruptHandler(int sig){
 
    int i;
    for(i = 0; i < 100 - primesCleared; i++){
        printf("fail\n\n");
    }  
    exit(0);
}
 
int main(void){
    
    int pid;
    
    /*pid = fork();
    primesCleared = 0;
    // Child process should kill execution after 14.9 s
    if(pid == 0){
        usleep(14900000);
        //printf("Killing!\n");
        kill(getppid(), SIGINT);
        return EXIT_SUCCESS;
    }
    else{
        signal(SIGINT, interruptHandler);*/
        nrOfFactors = 0;
        factors = (mpz_t*)malloc(sizeof(mpz_t)*150);
        initArray(factors, 150);
     
        mpz_t compositeNumber;
        mpz_init(compositeNumber);
        
        char inputBuffer[150];
     
        // Read number to be factored till there aren't any left
        while(fgets(inputBuffer, 150, stdin)){
            mpz_set_str(compositeNumber, inputBuffer, 10);  
            if(factorThis(compositeNumber)){
                //printFactors();
            }
            else{
                printf("fail\n\n");
            }
            primesCleared++;
            nrOfFactors = 0;
        }
     
        return EXIT_SUCCESS;
    
 
}

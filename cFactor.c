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
#define RHO 1 

mpz_t* factors;
int nrOfFactors;
mpz_t limit;
int primesCleared;
 
void addToFactors(mpz_t factor){
 
    mpz_set(factors[nrOfFactors], factor);
    nrOfFactors++;
 
}
 
void trialDivision(mpz_t compositeNumber){
 
    mpz_t test_prime;
    mpz_t upper_bound;
    mpz_t mod_test_prime;
 
    mpz_init(test_prime);
    mpz_init(upper_bound);
    mpz_init(mod_test_prime);
 
    mpz_set_ui(test_prime, 2);
    mpz_set_ui(upper_bound, 3587);
 
    while(mpz_cmp(test_prime, upper_bound) <= 0){
        //TODO square root something
        if(mpz_probab_prime_p(compositeNumber, 2)){
            return;
        }
 
        mpz_mod(mod_test_prime, compositeNumber, test_prime);
        if(mpz_cmp_ui(mod_test_prime, 0) == 0){
            addToFactors(test_prime);
            mpz_fdiv_q(compositeNumber, compositeNumber, test_prime);
        }
        else{
            mpz_nextprime(test_prime, test_prime);
        }
    }
 
}
 
void f_rho(mpz_t x, mpz_t n, unsigned int a){
 
    mpz_pow_ui(x, x, 2);
    mpz_add_ui(x, x, a);
    mpz_mod(x, x, n);
 
}
 
void pollardRhoStd(mpz_t compositeNumber, mpz_t res){
 
    unsigned int j = 0;
    mpz_t x;
    mpz_t y;
    mpz_t d;
    
    mpz_init(x);
    mpz_init(y);
    mpz_init(d);
 
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
 
    //205000 best
    while(j < 205000){
        f_rho(x, compositeNumber, 1);
        f_rho(y, compositeNumber, 1);
        f_rho(y, compositeNumber, 1);
    
        mpz_sub(d, x, y);
        mpz_abs(d, d);
        mpz_gcd(d, d, compositeNumber);
 
        j++;
        if(mpz_cmp_ui(d, 1) != 0){
            if(mpz_cmp(d, compositeNumber) != 0){
                mpz_set(res, d);
                return;
            }
            else{
                mpz_set_ui(res, 0);
                return;
            }
        }
        
    }
    
    mpz_set_ui(res, 0);
    return;
 
}

void min(mpz_t a, mpz_t b, mpz_t res){

    if(mpz_cmp(a, b) < 0)
        mpz_set(res, a);
    else
        mpz_set(res, b);

}

void pollardRhoBrent(mpz_t compositeNumber, mpz_t res){

    mpz_t x;
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
    // m = 10, fast and two erros m = 5 a little bit slower zero errors
    mpz_set_ui(m, 100);
    mpz_set_ui(q, 1);
    mpz_set_ui(ys, 1);
    mpz_set_ui(temp, 1);
    unsigned int i;

    int limitC = 0;
    int limit = 34500;
    while(mpz_cmp_ui(d, 1) == 0 && limitC < limit){
        mpz_set(x, y);
        
        for(i = 0; mpz_cmp_ui(r, i) > 0; i++){
            f_rho(y, compositeNumber, RHO);     
        }
        
        mpz_set_ui(k, 0);
       
        while((mpz_cmp(k, r) < 0) && (mpz_cmp_ui(d, 1) == 0) && (limitC < limit)){
            mpz_set(ys, y);
            mpz_sub(minVal, r, k);
            min(m, minVal, minVal);
            
            for(i = 0; mpz_cmp_ui(minVal, i) > 0; i++){
                f_rho(y, compositeNumber, RHO);
                mpz_sub(temp, x, y);
                mpz_abs(temp, temp);
                mpz_mul(q, temp, q);
                mpz_mod(q, q, compositeNumber);
            }
            
            mpz_gcd(d, q, compositeNumber);
            mpz_add(k, k, m);
                        limitC++;
        }
        
        mpz_mul_ui(r, r, 2);
    }
     
     if(limitC == limit){
        mpz_set_ui(res,0);
        return;
    }

    if(mpz_cmp(r, compositeNumber) == 0){
        limitC = 0;
        while(TRUE){
            f_rho(ys, compositeNumber, RHO);
            mpz_sub(temp, x, ys);
            mpz_abs(temp, temp);
            mpz_gcd(d, temp, compositeNumber);
            if(mpz_cmp_ui(d, 1) != 0 && limitC < 1000)
                break;
            limitC++;       
        }
    }
      
    if(mpz_cmp(r, compositeNumber) == 0 || limitC == 1000)
        mpz_set(res, 0);
    else
        mpz_set(res, d);
    
    return;
}

int factorThis(mpz_t compositeNumber){
    
    trialDivision(compositeNumber);
    mpz_t res;
    mpz_init(res);
    while(TRUE){
        if(mpz_probab_prime_p(compositeNumber, 2)){
            addToFactors(compositeNumber);
            return TRUE;
        }
 
        /*if(mpz_cmp(compositeNumber, limit) > 0){
            return FALSE;    
        }*/


        //pollardRhoStd(compositeNumber, res);
        pollardRhoBrent(compositeNumber, res);     

        if(mpz_cmp_ui(res, 0) == 0){
            return FALSE;    
        }
 
        if(mpz_probab_prime_p(res, 2)){
            addToFactors(res);    
        }
        else
            return FALSE;
 
        mpz_fdiv_q(compositeNumber, compositeNumber, res);
    }
    return FALSE;
}
 
void initArray(mpz_t* array, int size){
 
    int i;
    for(i = 0; i < size; i++){
        mpz_init(array[i]);
    }
 
}
 
void printFactors(){
    
    int i;
    for(i = 0; i < nrOfFactors; i++){
        gmp_printf("%Zd\n", factors[i]);
        /*if(!mpz_probab_prime_p(factors[i], 32)){
            printf("NO PRIME!\n");
        } */   
    }
    printf("\n");
}

void interruptHandler(int sig){
 
    //printf("Got killed, primes cleared so far: %d\n", primesCleared);
    int i;
    for(i = 0; i < 100 - primesCleared; i++){
        printf("fail\n\n");
    }  
    exit(0);
}
 
int main(void){
    
    int pid;
    
    pid = fork();
    primesCleared = 0;
    if(pid == 0){
        usleep(14900000);
        //printf("Killing!\n");
        kill(getppid(), SIGINT);
        return EXIT_SUCCESS;
    }
    else{
        //signal(SIGINT, interruptHandler);
        nrOfFactors = 0;
        factors = (mpz_t*)malloc(sizeof(mpz_t)*150);
        initArray(factors, 150);
        
        mpz_init(limit);
                            
                            //158456325028528675187087900672  97  bits best
                            //2535301200456458802993406410752 100 bits
                            //19807040628566084398385987584   93  bits
                            //2475880078570760549798248448    90  bits
                            //650000000000000000000000000    ~90  bits java
                            //77371252455336267181195264      85  bits
        mpz_set_str(limit, "650000000000000000000000000", 10);
     
        mpz_t compositeNumber;
        mpz_init(compositeNumber);
        
        char inputBuffer[150];
     
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
     
        //factorThis(test);
        //printFactors();
     
        return EXIT_SUCCESS;
    }
 
}

# 21/11/2018
# k-hit
# + Each mutation has an ID (ID) and a growth advantage (s) as a tuple
# ``[(ID1, s1)]``
# + Each clone contains a certain number of mutations, so could be defined as a list of the mutation, fitness effect tuples
# ``[(ID1, s1),(ID2, s2)....]``
# + Each clone also needs an associated size (n) (e.g. 50 cells) -> ``[[(ID1, s1),(ID2, s2)....], n]``
# + Total population is a list of the different clones (with a different set of mutations within it): ``[[[(ID1, s1),(ID2, s2)....], n], [[(ID1, s1),(ID2, s2)....], n]]]``
#this turns on mean fitness update

#import libraries
library(rlist)
library(distr)
library(ggplot2)
library(scales)
library(gridExtra)
library(beepr)


start_time <- Sys.time()

#par(mfrow=c(3,1))

# Bio_system <-list(list(list(list(1,0)),10^5));#clone_info<- (mutation, clone size)
#starting from one single clone with three mutations
#Bio_system <-list(list(list(list(1,0.01),list(2,0.02),list(3,0.0015)),100),list(list(list(4,0.1)),100000));#clone_info<- (mutation, clone size)
    #starting from one single clone with three mutations

#Define function 'Divide'
Divide<- function(clone_info, dt, B0, D0,mean_fitness){
  #dt = a small interval of time, 
  #D0 = symmetric rate to differentiated
  #B0 = symmetric rate to self renewal
  #B  = modified birth rate including fitness advantage
  
      cloneID=clone_info[[1]] # [[ID1, s1],[ID2,s2]...] is the ID of the clone
      clone_size=clone_info[[2]][1] # n
#     print('clone_size =',clone_size)
      
      Total_fitness_effects=0;
      for (i in 1:length(cloneID) )
        {Total_fitness_effects=Total_fitness_effects+cloneID[[i]][[2]]}
      #fitness of the clone is the sum of fitness conferred by all mutations it contains
      Total_fitness_effects=Total_fitness_effects- mean_fitness                                     

B = B0+Total_fitness_effects #B = birth rate, which is 1 + the overall fitness effect of the mutations within a particular clone
if (clone_size!=0){
  if(B>0){
number_of_births <- rpois(1, clone_size*B*dt)#pulls a random number from the Poisson distribution with/
                                                #mean of clone_size x birth rate x interval of time

          }else{number_of_births=0}            #to make sure the mean is positive        

number_of_deaths <- rpois(1, clone_size*D0*dt);
#     print('number of deaths = ', number_of_deaths)
}
        else{number_of_births=0;
              number_of_deaths=0; }          #makes sure an extinct clone stays extinct


new_clone_size = clone_size + number_of_births - number_of_deaths;
if(new_clone_size<0){new_clone_size=0}
      #print(paste("the new clone size is ",new_clone_size))
return (new_clone_size)
}


#Define function 'mutation_fitness' 
mutation_fitness <- function(DFE,ratio_neutral_beneficial,benefit) {
  # library(distr)
  # DFE    <- function(x) (2/pi) * (1/(exp(x)+exp(-x)))  # Define the continuous probability density function DFE
  # dist_DFE <-distr::AbscontDistribution(d=DFE)  # signature for a dist with DFE 
  # fitness <- distr::r(dist_DFE) 
  # fitness_effect=fitness(1)
  # if(abs(fitness_effect)>0.1){fitness_effect=0.1} #to restrict the range of fitness_effect's
  # 
  random_number<- runif(1, 0, 1)
  threshold= 1/(1+ratio_neutral_beneficial)
  if(random_number<threshold){
    
    s=benefit
    
  }else{s=0}
  return(s)                                  
}


#Define function 'mutate'
mutate<-function(clone_info, dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit){ #u=mutation rate
  mutations=clone_info[[1]] #[[ID1, s1],[ID2,s2]]
  clone_size=clone_info[[2]]
  
  # if(clone_size*u*dt>1){                                       #do not mutate clones with size<1/U
  
  
  number_of_mutations <-rpois(1, clone_size*u*dt) ;
     #print(paste("mean of poisson =", clone_size*u*dt))
      # print(paste("number of mutations =", number_of_mutations))

list_of_new_clones <- list()
new_clone_info <-list()
if(number_of_mutations!=0){
for (i in 1:number_of_mutations){
  last_mutation_ID = last_mutation_ID+1
        #print(paste("last mutation ID =", last_mutation_ID))
  new_fitness_effect = mutation_fitness(1,ratio_neutral_beneficial,benefit)
#         print("new fitness effect =", new_fitness_effect)
  new_mutation<-list(last_mutation_ID, new_fitness_effect)
  new_clone_info <-list(list.append(mutations,new_mutation),1)    #[[[ID1,s1],[ID2,s2],....],n]
#         print("new clone =", new_clone)
  list_of_new_clones<-list.append(list_of_new_clones,new_clone_info)
#         print("list of new clones =", list_of_new_clones)
}
}   
   
   
   
  # }else{
  #   
  #   list_of_new_clones<-list()
  #   last_mutation_ID = last_mutation_ID
  #   number_of_mutations=0
  #   
  # }
  
result<-list(list_of_new_clones, last_mutation_ID, number_of_mutations)     
return (result)                                      #result[1] returns the list of new clones
}
#what happens when number_of_mutations=0?

#define function 'update_mean_fitness'
update_mean_fitness<- function(Bio_system)
{        
        
        number_of_clones= length(Bio_system)
        Total_fitness_effects=0;
        total_pop=0;
        weighted_fitness=0;
         for (i in 1:number_of_clones){                #i is the tag of the clone
           number_of_mutations=length(Bio_system[[i]][[1]])  #number of mutations in that clone
           Total_fitness_effects=0
           for (j in 1:number_of_mutations) {
             Total_fitness_effects= Total_fitness_effects+Bio_system[[i]][[1]][[j]][[2]]
           }                                       #sum up the fitness conferred by each mutation in a clone
           
           weighted_fitness=weighted_fitness+Bio_system[[i]][[2]]*Total_fitness_effects
           total_pop=total_pop+Bio_system[[i]][[2]];                 
             
         }
        mean_fitness=weighted_fitness/total_pop;

return (list(mean_fitness, number_of_clones,total_pop))       
}


#Main function


Number_of_Participants=600000
Population_size=10^5                                 #is the starting population size of clone
hitchhiker_cutoff=0.5;
benefit=0.1;
delta_time=0.1;#in units of generation
number_of_generations=150;#number of generations
lifespan=number_of_generations/delta_time; #number of runs

N=1
VAFlim=0.1
smallest_detectable_size=2*VAFlim*Population_size



#initialise the list in N
Collective_Number_of_detectable_single_mutants<-list()
Collective_Number_of_detectable_double_mutants<-list()
for (t in 1:lifespan){
  Collective_Number_of_detectable_single_mutants[[N]]=0
  Collective_Number_of_detectable_double_mutants[[N]]=0
}

Clonal_histories_per_person<-list()
Entire_experiment<-list()

for (N in 1:Number_of_Participants) {
  
Bio_system <-list(list(list(list(1,0)),Population_size));
#population_limit=10^15;
t=0;
last_mutation_ID=1;
#*delta_time   is introduced in dt
dt=1*delta_time;
B0=1;
D0=1;
u=10^-5;                        #mutation rate (neutral and beneficial)
ratio_neutral_beneficial=10^-5;           #R=1 for beneficial mutation
mean_fitness=0;
Clonal_histories<-list() 

Stored_mean_fitness<-list()
Stored_number_of_clones<-list()
Bio_system_record<-list()
mean_fitness_record<-list()
Number_of_nonextinct_clones_record<-list()

#initialise the list in t
Number_of_detectable_single_mutants<-list()
Number_of_detectable_double_mutants<-list()
for (t in 1:lifespan){
Number_of_detectable_single_mutants[[t]]=0
Number_of_detectable_double_mutants[[t]]=0
}

for (t in 1:lifespan) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  for (m in 1:length(Bio_system)){
    
    if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
      if(length(Bio_system[[m]][[1]])<3){ 
      
      
      result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
      list_of_new_clones[[m]]<- result[[m]][[1]]
      last_mutation_ID <- result[[m]][[2]]
      number_of_mutations[[m]]<- result[[m]][[3]]
      #only call the function 'mutate' once for each clone
      }
    }
  }
  for (m in 1:length(Bio_system)){
    if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
      if(length(Bio_system[[m]][[1]])<3){ 
      
      if (number_of_mutations[[m]]>0){
        for (i in 1:number_of_mutations[[m]]) {
          
          if(length(result[[m]][[1]])!=0){
            Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
          }
          
          
        } 
      }
    }
    }
  }
  #print("Finished mutating")
  #divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  for (i in 1:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }
  #print("Finished dividing")
  
  #update mean fitness
  result2<-update_mean_fitness(Bio_system)
  mean_fitness<-result2[[1]]
  Stored_mean_fitness<-list.append(Stored_mean_fitness,mean_fitness)
  
  total_pop<-result2[[2]]
  
  
  #restrict population size : no need
  
  
  # Clonal_histories <- list.append(Clonal_histories,Clonal_sizes) 
  # # is a list of clonal sizes at different times
  # 
  # 
  # 
  # number_of_nonextinct_clones=0
  # for (m in 1:length(Bio_system)) {
  #   
  # if(Bio_system[[m]][[2]]!=0){
  # number_of_nonextinct_clones=number_of_nonextinct_clones+1
  #                            }
  # 
  # 
  # }
  # Stored_number_of_clones<-list.append(Stored_number_of_clones,number_of_nonextinct_clones)
  # Bio_system_record[[t]]<-Bio_system
  # #this stores the Bio_system at each time point into a separate entry for a single patient



#check not dead check length of clonID is 2
  nu=0
  mu=0
for (n in 1:length(Bio_system)) {
  if (n!=1){
  CloneID<- Bio_system[[n]][[1]]
  Size<- Bio_system[[n]][[2]]
  if (Size>smallest_detectable_size){
    if(length(CloneID)==2){nu=nu+1}
    
    if(length(CloneID)==3){mu=mu+1}
  }
}
}


  
  t=t+dt;
  #print(paste("this is generation ",t))
  # this is logic 1 verus 0: Yes or no
  if (nu>=2){Number_of_detectable_single_mutants[[t]]=1 }
  if (mu>=1){Number_of_detectable_double_mutants[[t]]=1 }
  
}

 




Collective_Number_of_detectable_double_mutants[[N]]= Number_of_detectable_double_mutants
Collective_Number_of_detectable_single_mutants[[N]]= Number_of_detectable_single_mutants



                   print(paste("This is patient number ",N, " has number of mutations", last_mutation_ID ))
                   
                   
                   
                   
}                                                  #end of number of participant loop


end_time <- Sys.time()

code_ran_for<- end_time - start_time

code_ran_for<-round(code_ran_for,1)

rate_detecting_single_mutant<-list()
rate_detecting_double_mutant<-list()
for (t in 1:lifespan){
  rate_detecting_single_mutant[[t]]=0
  rate_detecting_double_mutant[[t]]=0
        for (N in 1:Number_of_Participants){
          rate_detecting_single_mutant[[t]]=rate_detecting_single_mutant[[t]]+Collective_Number_of_detectable_single_mutants[[N]][[t]]/Number_of_Participants
          rate_detecting_double_mutant[[t]]=rate_detecting_double_mutant[[t]]+Collective_Number_of_detectable_double_mutants[[N]][[t]]/Number_of_Participants
          
        }

}


write.csv(unlist(rate_detecting_single_mutant),file='single_mutant_detection_rate_VAFlim_0.1_s_0.1_u_minus5_600000ppl.csv')
write.csv(unlist(rate_detecting_double_mutant),file='double_mutant_detection_rate_VAFlim_0.1_s_0.1_u_minus5_600000ppl.csv')

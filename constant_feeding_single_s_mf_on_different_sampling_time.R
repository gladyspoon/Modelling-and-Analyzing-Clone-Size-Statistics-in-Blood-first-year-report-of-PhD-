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


#Define function 'Divide_develop'
Divide_develop<- function(clone_info, dt, B0, D0,increase,mean_fitness){
  
  cloneID=clone_info[[1]] # [[ID1, s1],[ID2,s2]...] is the ID of the clone
  clone_size=clone_info[[2]][1] # n
  # #     print('clone_size =',clone_size)
  # Total_fitness_effects=0
  # 
  # for (i in 1:length(cloneID) )
  # {Total_fitness_effects=Total_fitness_effects+cloneID[[i]][[2]]}
  # #fitness of the clone is the sum of fitness conferred by all mutations it contains
  # Total_fitness_effects=Total_fitness_effects- mean_fitness                                     
  # 
  # # #B=r+Total_fitness_effects
  # if (clone_size!=0){
  # 
  #     number_of_births <- ceiling(clone_size*r*dt)#rpois(1, clone_size*B*dt)#pulls a random number from the Poisson distribution with/
  #     #mean of clone_size x birth rate x interval of time
  #     
  #   
  # }
  # else{number_of_births=0;
  #  }          #makes sure an extinct clone stays extinct
  # 
  # 
  # #new_clone_size = clone_size + number_of_births;
  new_clone_size = clone_size+increase*dt
  if(new_clone_size<0){new_clone_size=0}
  
  return (N)
}


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
      Total_fitness_effects=Total_fitness_effects -mean_fitness                                     

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
  
  
  # random_number<- runif(1, 0, 1)
  # threshold= 1/(1+ratio_neutral_beneficial)
  # if(random_number<threshold){
  #   
  #   s=benefit
  #   
  # }else{s=0}
  
  return(benefit)                                  
}


#Define function 'mutate'
mutate<-function(clone_info, dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit){ #u=mutation rate
  mutations=clone_info[[1]] #[[ID1, s1],[ID2,s2]]
  clone_size=clone_info[[2]]
  
  # if(clone_size*u*dt>1){                                       #do not mutate clones with size<1/U
  
  
  number_of_mutations <-rpois(1, clone_size*u*dt) ;
     #print(paste("mean of poisson =", clone_size*u*dt))
      # print(paste("number of mutations =", number_of_mutations))

list_of_new_clones=list()
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
Number_of_Participants=50000
N=10^5
#Ideal_Population_size=10^7                              
delta_time=0.1;#in units of generation
#r=log(Ideal_Population_size)/stop6;
increase=16144.04 #increase in feeding population numbers per generation
 #growth rate of feeding clone
alpha=0;
r=0.15;
benefit=0.1; #is the fitness of the first mutant i.e. r-s_del
#this formalism means double mutant have a fitness=2*benefit which is not 2 times the fitness conferred by the mutation
stop1=40;
stop2=60;
stop3=80;
stop4=90;
stop5=100;
stop6=120;



lifespan1=stop1/delta_time;
lifespan2=stop2/delta_time;
lifespan3=stop3/delta_time;
lifespan4=stop4/delta_time;
lifespan5=stop5/delta_time;
lifespan6=stop6/delta_time;
t_stop_mutating=lifespan6


Total_mutation_frequency_all_deleterious_1<-list()
          Total_mutation_frequency_all_deleterious_2<-list()
          Total_mutation_frequency_all_deleterious_3<-list()
          Total_mutation_frequency_all_deleterious_4<-list()
          Total_mutation_frequency_all_deleterious_5<-list()
          Total_mutation_frequency_all_deleterious_6<-list()
Total_mutation_frequency_single_deleterious_1<-list()
          Total_mutation_frequency_single_deleterious_2<-list()
          Total_mutation_frequency_single_deleterious_3<-list()
          Total_mutation_frequency_single_deleterious_4<-list()
          Total_mutation_frequency_single_deleterious_5<-list()
          Total_mutation_frequency_single_deleterious_6<-list()
Clonal_histories_per_person<-list()
Entire_experiment<-list()

for (N in 1:Number_of_Participants) {
  
Bio_system <-list(list(list(list(1,0)),N));
#population_limit=10^15;
t=0;
last_mutation_ID=1;
#*delta_time   is introduced in dt
dt=1*delta_time;
B0=1;
D0=1;
u=10^-5;                        #deleterious mutation rate 
ratio_neutral_beneficial=10/1;           #is irrelevant in this code
mean_fitness=0;

Clonal_histories<-list() 
Stored_mean_fitness<-list()
Stored_number_of_clones<-list()
Bio_system_record<-list()
mean_fitness_record<-list()
Number_of_nonextinct_clones_record<-list()


for (t in 1:lifespan1) {


  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
  #print("Finished mutating")
#divide: propagate growth of each clone

 
          Clonal_sizes<-list()
          
          #feeding clone
          Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
           
            
                                                                 #consider relative fitness
          Clonal_sizes[[1]]=Bio_system[[1]][[2]];
          
          #for all the rest mutant clones
          if(length(Bio_system)!=1){
          for (i in 2:length(Bio_system)) {
            
            if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
            
              if(TRUE){
            Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
                   }    #else{Bio_system[[i]][[2]]=population_limit}
            
            }                                                       #consider relative fitness
            Clonal_sizes[[i]]=Bio_system[[i]][[2]];
          
          }}
  #print("Finished dividing")
          
#update mean fitness
          result2<-update_mean_fitness(Bio_system)
          mean_fitness<-result2[[1]]
          Stored_mean_fitness<-list.append(Stored_mean_fitness,mean_fitness)
          
          total_pop<-result2[[2]]



         
          
          #  Clonal_histories <- list.append(Clonal_histories,Clonal_sizes) 
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
           #this stores the Bio_system at each time point into a separate entry for a single patient
          #print(paste("this is ",t*delta_time,"-th generation"))
          t=t+dt;
          #print(paste("this is generation ",t))
         
}



#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     


#######################################START OF PART FOR CLONAL HISTORY#################################
#to make all vectors of Clonal_histories the same length
# T<- lifespan1-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 

# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
# 
#            for (t in 1:lifespan1){
# 
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
# 
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
# #stores clonal histories of each person into a separate entry
# 
# 
# 
# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_1<-list()                              #all deleterious
mutation_frequency_single_deleterious_1<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_1[[m]]=0;
                              mutation_frequency_single_deleterious_1[[m]]=0;
                              
                              }    #initialization

for (m in 2:last_mutation_ID) {
 
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
   #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_1[[m]]=mutation_frequency_all_deleterious_1[[m]]+Size
          if(total_fitness==benefit){
        mutation_frequency_single_deleterious_1[[m]]=mutation_frequency_single_deleterious_1[[m]]+Size
          }
             
        
        }
            
        }       
      
      
      
    }
    
  }
 
 
}




Total_mutation_frequency_all_deleterious_1<-list.append(Total_mutation_frequency_all_deleterious_1,mutation_frequency_all_deleterious_1)
Total_mutation_frequency_single_deleterious_1<-list.append(Total_mutation_frequency_single_deleterious_1,mutation_frequency_single_deleterious_1)



# resume

for (t in lifespan1:lifespan2) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
  #print("Finished mutating")
  
#divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  #feeding clone
  Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
  
  
  #consider relative fitness
  Clonal_sizes[[1]]=Bio_system[[1]][[2]];
  
  #for all the rest mutant clones
  if(length(Bio_system)!=1){
  for (i in 2:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }}
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
  #print(paste("this is ",t*delta_time,"-th generation"))
  t=t+dt;
  #print(paste("this is generation ",t))
  
}


#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     



#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_2<-list()                              #all deleterious
mutation_frequency_single_deleterious_2<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_2[[m]]=0;
mutation_frequency_single_deleterious_2[[m]]=0;

}    #initialization

for (m in 2:last_mutation_ID) {
  
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
    #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_2[[m]]=mutation_frequency_all_deleterious_2[[m]]+Size
          if(total_fitness==benefit){
            mutation_frequency_single_deleterious_2[[m]]=mutation_frequency_single_deleterious_2[[m]]+Size
          }
          
          
        }
        
      }       
      
      
      
    }
    
  }
  
  
}




Total_mutation_frequency_all_deleterious_2<-list.append(Total_mutation_frequency_all_deleterious_2,mutation_frequency_all_deleterious_2)
Total_mutation_frequency_single_deleterious_2<-list.append(Total_mutation_frequency_single_deleterious_2,mutation_frequency_single_deleterious_2)







#resume once again
# resume

for (t in lifespan2:lifespan3) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
#divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  #feeding clone
  Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
  
  
  #consider relative fitness
  Clonal_sizes[[1]]=Bio_system[[1]][[2]];
  
  #for all the rest mutant clones
  if(length(Bio_system)!=1){
  for (i in 2:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }}
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
  #print(paste("this is ",t*delta_time,"-th generation"))
  t=t+dt;
  #print(paste("this is generation ",t))
  
}


#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     



#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_3<-list()                              #all deleterious
mutation_frequency_single_deleterious_3<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_3[[m]]=0;
mutation_frequency_single_deleterious_3[[m]]=0;

}    #initialization

for (m in 2:last_mutation_ID) {
  
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
    #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_3[[m]]=mutation_frequency_all_deleterious_3[[m]]+Size
          if(total_fitness==benefit){
            mutation_frequency_single_deleterious_3[[m]]=mutation_frequency_single_deleterious_3[[m]]+Size
          }
          
          
        }
        
      }       
      
      
      
    }
    
  }
  
  
}




Total_mutation_frequency_all_deleterious_3<-list.append(Total_mutation_frequency_all_deleterious_3,mutation_frequency_all_deleterious_3)
Total_mutation_frequency_single_deleterious_3<-list.append(Total_mutation_frequency_single_deleterious_3,mutation_frequency_single_deleterious_3)





for (t in lifespan3:lifespan4) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
  #print("Finished mutating")
  
#divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  #feeding clone
  Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
  
  
  #consider relative fitness
  Clonal_sizes[[1]]=Bio_system[[1]][[2]];
  
  #for all the rest mutant clones
  if(length(Bio_system)!=1){
  for (i in 2:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }}
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
  #print(paste("this is ",t*delta_time,"-th generation"))
  t=t+dt;
  #print(paste("this is generation ",t))
  
}


#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     



#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_4<-list()                              #all deleterious
mutation_frequency_single_deleterious_4<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_4[[m]]=0;
mutation_frequency_single_deleterious_4[[m]]=0;

}    #initialization

for (m in 2:last_mutation_ID) {
  
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
    #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_4[[m]]=mutation_frequency_all_deleterious_4[[m]]+Size
          if(total_fitness==benefit){
            mutation_frequency_single_deleterious_4[[m]]=mutation_frequency_single_deleterious_4[[m]]+Size
          }
          
          
        }
        
      }       
      
      
      
    }
    
  }
  
  
}




Total_mutation_frequency_all_deleterious_4<-list.append(Total_mutation_frequency_all_deleterious_4,mutation_frequency_all_deleterious_4)
Total_mutation_frequency_single_deleterious_4<-list.append(Total_mutation_frequency_single_deleterious_4,mutation_frequency_single_deleterious_4)





for (t in lifespan4:lifespan5) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
  #print("Finished mutating")

  #divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  #feeding clone
  Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
  
  
  #consider relative fitness
  Clonal_sizes[[1]]=Bio_system[[1]][[2]];
  
  #for all the rest mutant clones
  if(length(Bio_system)!=1){
  for (i in 2:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }}
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
  #print(paste("this is ",t*delta_time,"-th generation"))
  t=t+dt;
  #print(paste("this is generation ",t))
  
}

#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     




#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_5<-list()                              #all deleterious
mutation_frequency_single_deleterious_5<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_5[[m]]=0;
mutation_frequency_single_deleterious_5[[m]]=0;

}    #initialization

for (m in 2:last_mutation_ID) {
  
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
    #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_5[[m]]=mutation_frequency_all_deleterious_5[[m]]+Size
          if(total_fitness==benefit){
            mutation_frequency_single_deleterious_5[[m]]=mutation_frequency_single_deleterious_5[[m]]+Size
          }
          
          
        }
        
      }       
      
      
      
    }
    
  }
  
  
}




Total_mutation_frequency_all_deleterious_5<-list.append(Total_mutation_frequency_all_deleterious_5,mutation_frequency_all_deleterious_5)
Total_mutation_frequency_single_deleterious_5<-list.append(Total_mutation_frequency_single_deleterious_5,mutation_frequency_single_deleterious_5)



print(paste("entering generation ", lifespan5*delta_time))

for (t in lifespan5:lifespan6) {
  
  
  #mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()
  
  #for (m in 1:length(Bio_system)){
  m=1
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
    list_of_new_clones[[m]]<- result[[m]][[1]]
    last_mutation_ID <- result[[m]][[2]]
    number_of_mutations[[m]]<- result[[m]][[3]]
    #only call the function 'mutate' once for each clone
    
  }
  #}
  
  if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
    
    if (number_of_mutations[[m]]>0){
      for (i in 1:number_of_mutations[[m]]) {
        
        if(length(result[[m]][[1]])!=0){
          Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
        }
        
        
      } 
    }
    
  }
  #print("Finished mutating")
  
  #divide: propagate growth of each clone
  
  
  Clonal_sizes<-list()
  
  #feeding clone
  Bio_system[[1]][[2]]=Divide_develop(Bio_system[[1]], dt, B0, D0,increase,mean_fitness)[[1]]   #update feeding clone size
  
  
  #consider relative fitness
  Clonal_sizes[[1]]=Bio_system[[1]][[2]];
  
  #for all the rest mutant clones
  if(length(Bio_system)!=1){
  for (i in 2:length(Bio_system)) {
    
    if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
      
      if(TRUE){
        Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
      }    #else{Bio_system[[i]][[2]]=population_limit}
      
    }                                                       #consider relative fitness
    Clonal_sizes[[i]]=Bio_system[[i]][[2]];
    
  }}
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
  #print(paste("this is ",t*delta_time,"-th generation"))
  t=t+dt;
  #print(paste("this is generation ",t))
  
}

#End_Population_Size=1
#for(m in 1:length(Bio_system)){End_Population_Size=End_Population_Size+Bio_system[[m]][[2]]}     




#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_all_deleterious_6<-list()                              #all deleterious
mutation_frequency_single_deleterious_6<-list()


for (m in 1:last_mutation_ID){mutation_frequency_all_deleterious_6[[m]]=0;
mutation_frequency_single_deleterious_6[[m]]=0;

}    #initialization

for (m in 2:last_mutation_ID) {
  
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    
    #make the lists for mutation frequencies
    
    for (k in 1:length(CloneID)) {
      
      if(CloneID[[k]][[1]]==m){
        
        if(CloneID[[k]][[2]]==benefit){
          mutation_frequency_all_deleterious_6[[m]]=mutation_frequency_all_deleterious_6[[m]]+Size
          if(total_fitness==benefit){
            mutation_frequency_single_deleterious_6[[m]]=mutation_frequency_single_deleterious_6[[m]]+Size
          }
          
          
        }
        
      }       
      
      
      
    }
    
  }
  
  
}




Total_mutation_frequency_all_deleterious_6<-list.append(Total_mutation_frequency_all_deleterious_6,mutation_frequency_all_deleterious_6)
Total_mutation_frequency_single_deleterious_6<-list.append(Total_mutation_frequency_single_deleterious_6,mutation_frequency_single_deleterious_6)



                   print(paste("This is patient number ",N, " has number of mutations", last_mutation_ID ))
                   
                   
                   
                   
}                                                 #end of number of participant loop



end_time <- Sys.time()

code_ran_for<- end_time - start_time

code_ran_for<-round(code_ran_for,1)
#write.csv(unlist(Total_mutation_frequency_single_deleterious_2), file= "linear_feeding_increase_16144_s2_0.1_5000ppl_age_110.csv")
#write.csv(unlist(Total_mutation_frequency_single_deleterious_4), file= "alpha_0_5000ppl_Age90_u_10minus6")
#write.csv(unlist(Total_mutation_frequency_single_deleterious_6), file= "alpha_0_5000ppl_Age110_u_10minus6")



######################################################################################################################
###################################     PLOT HISTORY OF A SINGLE PATIENT     ###################################


# 
# M=1     #specify which patient's history
# 
# #plot the clonal histories
# 
# History<-Clonal_histories_per_person[[M]]
# timescale<-c(1:lifespan1)
# matplot(timescale, log(History[,1:length(Bio_system)]),col=1:length(Bio_system),lty=2, lwd=1,
#         type = "l",
#         xlab="number of generations",
#         ylab="Clone Size",
#         xlim=c(0,lifespan1)
#        #ylim=c(0,110000)
# 
# )
# #title(main = "Trajectories of mutants, 1000 generations, 8000 runs")
# 
# 
# 
# #plot mean fitness with time
# 
# Record_mean_fitness<-unlist(mean_fitness_record[[M]])
# 
# plot(timescale,Record_mean_fitness,xlab="number of generations",
#      ylab="Mean fitness",lty=2)
# 
# 
# 
# #plot mean fitness with time
# 
# Number_of_nonextinct_clones<-unlist(Number_of_nonextinct_clones_record[[M]])
# 
# plot(timescale,Number_of_nonextinct_clones,xlab="number of generations",
#      ylab="Number of nonextinct Clones",lty=2)
# 
# 
# 
# 
# grid.arrange(plot1,plot2,plot3, nrow = 3)
# 



##################################    OVERVIEW of ALL PATIENTS    ###################################################




# #Mutation frequency distribution
# Sum<-list.append(Total_mutation_frequency_neutral,Total_mutation_frequency_beneficial)
# Sum<-unlist(Sum)
# total<-Sum[Sum!=0]
# total<-total/Population_size                                            #normalize total mutation frequency
# logtotal<-log(total)
# 
# Total_mutation_frequency_beneficial<-unlist(Total_mutation_frequency_beneficial)
# Total_mutation_frequency_beneficial_first_mutant<-unlist(Total_mutation_frequency_beneficial_first_mutant)
# Total_mutation_frequency_beneficial_double_mutant<-unlist(Total_mutation_frequency_beneficial_double_mutant)
# Total_mutation_frequency_neutral<-unlist(Total_mutation_frequency_neutral)
# Total_mutation_frequency_neutral_not_hitchhiker<-unlist(Total_mutation_frequency_neutral_not_hitchhiker)
# Total_mutation_frequency_neutral_hitchhiker<-unlist(Total_mutation_frequency_neutral_hitchhiker)
# Total_mutation_frequency_neutral_beneficial_first<-unlist(Total_mutation_frequency_neutral_beneficial_first)
# Total_mutation_frequency_neutral_beneficial_later<-unlist(Total_mutation_frequency_neutral_beneficial_later)
# 
# p<-Total_mutation_frequency_beneficial[Total_mutation_frequency_beneficial!=0]
# p<-p/Population_size                                            #normalize beneficial mutation frequency
# logp<-log(p)
# 
# p2<-Total_mutation_frequency_beneficial_first_mutant[Total_mutation_frequency_beneficial_first_mutant!=0]
# p2<-p2/Population_size                         #normalize beneficial mutation frequency for those in first mutants
# logp2<-log(p2)
# 
# p3<-Total_mutation_frequency_beneficial_double_mutant[Total_mutation_frequency_beneficial_double_mutant!=0]
# p3<-p3/Population_size                         
# logp3<-log(p3)
# 
# 
# q<-Total_mutation_frequency_neutral[Total_mutation_frequency_neutral!=0]
# q<-q/Population_size                                            #normalize neutral mutation frequency
# logq<-log(q)
# 
# q2<-Total_mutation_frequency_neutral_not_hitchhiker[Total_mutation_frequency_neutral_not_hitchhiker!=0]
# q2<-q2/Population_size                        #normalize neutral mutation frequency for those that are not hitchhikers
# logq2<-log(q2)
# 
# q3<-Total_mutation_frequency_neutral_beneficial_first[Total_mutation_frequency_neutral_beneficial_first!=0]
# q3<-q3/Population_size                        #normalize neutral mutation frequency for type "beneficial first"
# logq3<-log(q3)
# 
# q4<-Total_mutation_frequency_neutral_beneficial_later[Total_mutation_frequency_neutral_beneficial_later!=0]
# q4<-q4/Population_size                        #normalize neutral mutation frequency for type "beneficial later"
# logq4<-log(q4)
# 
# q5<-Total_mutation_frequency_neutral_hitchhiker[Total_mutation_frequency_neutral_hitchhiker!=0]
# q5<-q5/Population_size                        #normalize neutral mutation frequency for all hitchhikers
# logq5<-log(q5)
# 
# 
# 
# 
# 
# BINWIDTH<-0.1
# lower_xlim=log(20/Population_size)                    #only show mutation frequency large enough to be continuous in bin sizes
# 
# u_neu=u*ratio_neutral_beneficial/(ratio_neutral_beneficial+1);
# u_ben=u/(ratio_neutral_beneficial+1);
# n_beneficial=(exp(benefit*lifespan)-1)/benefit
# n_neutral=lifespan
# n_hitchhiker=(exp(benefit*lifespan))/benefit
# 
# plot1<- qplot(logp, geom="histogram",
#       #log = "x",
#       
#       binwidth=(function(x) BINWIDTH),
#       xlab = "Variant Allele Frequency (x100 %)",
#       main = paste("All Beneficial mutations"
#         # #" s =", benefit,
#         #            "\n U =", u/(1+ratio_neutral_beneficial),  "          lifespan =",lifespan,"          population size =",
#         #            Population_size,"\n Number of Participants =",
#         #            Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#                    )) + scale_y_continuous(trans = log_trans(), 
#                                                                  breaks = trans_breaks("log", function(x) exp(x)),
#                                                                  labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity", 
#                                                                  breaks = c(-9,-7,-5,-3,-1,1),
#                                                                  labels =  trans_format("identity", math_format(e^.x)),
#                                                                  limits =  c(lower_xlim,1)
#                                                                  )+stat_function( 
#                                                                    fun = function(x){
#                                                                      #1/x
#                                                                      
#                                                                      BINWIDTH*Number_of_Participants*Population_size*u_ben*exp(-exp(x+log(Population_size))/n_beneficial)
#                                                                    },
#                                                                    xlim =  c(lower_xlim,-1.2)
#                                                                    #, 
#                                                                    #args = c(mean = mean, sd = sd, n = n, bw = binwidth)
#                                                                  )
# 
# plot2<- qplot(logp2, geom="histogram",
#               #log = "x",
#               
#               binwidth=(function(x) BINWIDTH),
#               xlab = "Variant Allele Frequency (x100 %)",
#               main = paste("First beneficial mutants"
#                            # #" s =", benefit,
#                            # "\n U =", u/(1+ratio_neutral_beneficial),  "          lifespan =",lifespan,"          population size =",
#                            # Population_size,"\n Number of Participants =",
#                            # Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#                            )) + scale_y_continuous(trans = log_trans(), 
#                                                                                                                         breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                                                         labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity", 
#                                                                                                                                                                                               breaks = c(-9,-7,-5,-3,-1,1),
#                                                                                                                                                                                               labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                                               limits =  c(lower_xlim,1)
#                                                                                                                         )
# 
# 
# 
# plot3<- qplot(logq, geom="histogram",
#       #log = "x",
#       
#       binwidth=(function(x) BINWIDTH),
#       xlab = "Variant Allele Frequency (x100 %)",
#       main = paste("All Neutral mutations"
#         # #" s =", benefit,
#         # "\n U =", u*(1-1/(1+ratio_neutral_beneficial)),  "          lifespan =",lifespan,"          population size =",
#         # Population_size,"\n Number of Participants =",
#         # Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#         )) + scale_y_continuous(trans = log_trans(), 
#                                                                                               breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                               labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity",
#                                                                                                                                                                     breaks = c(-9,-7,-5,-3,-1,1),
#                                                                                                                                                                     labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                     limits =  c(lower_xlim,1)
#                                                                                               )+stat_function( 
#                                                                                                 fun = function(x){
#                                                                                                   #1/x
#                                                                                                   #Population_size*u*exp(-exp(x)/n_beneficial)
#                                                                                                   BINWIDTH*Number_of_Participants*Population_size*exp(-exp(x+log(Population_size))/n_neutral)*u_neu
#                                                                                                 },
#                                                                                                 xlim =  c(lower_xlim,-6.5)
#                                                                                                 #, 
#                                                                                                 #args = c(mean = mean, sd = sd, n = n, bw = binwidth)
#                                                                                               )
# 
# 
# plot4<- qplot(logq2, geom="histogram",
#               #log = "x",
#               
#               binwidth=(function(x) BINWIDTH),
#               xlab = "Variant Allele Frequency (x100 %)",
#               main = paste("Non-hitchhikers"
#                            # #" s =", benefit,
#                            # "\n U =", u*(1-1/(1+ratio_neutral_beneficial)),  "          lifespan =",lifespan,"          population size =",
#                            # Population_size,"\n Number of Participants =",
#                            # Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#                            )) + scale_y_continuous(trans = log_trans(), 
#                                                                                                                         breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                                                         labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity",
#                                                                                                                                                                                               breaks = c(-9,-7,-5,-3,-1,1),
#                                                                                                                                                                                               labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                                               limits =  c(lower_xlim,1)
#                                                                                                                         )
# plot5<- qplot(logq3, geom="histogram",
#               #log = "x",
#               
#               binwidth=(function(x) BINWIDTH),
#               xlab = "Variant Allele Frequency (x100 %)",
#               main = paste("Hitchhikers with beneficial mutation coming in earlier"
#                            # #" s =", benefit,
#                            # "\n U =", u*(1-1/(1+ratio_neutral_beneficial)),  "          lifespan =",lifespan,"          population size =",
#                            # Population_size,"\n Number of Participants =",
#                            #Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#                            )) + scale_y_continuous(trans = log_trans(),
#                                                                                                                         breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                                                         labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity",
#                                                                                                                                                                                               breaks = c(-24,-20,-16,-12),
#                                                                                                                                                                                               labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                                               limits =  c(lower_xlim,1)
#                                                                                                                         )
# plot6<- qplot(logq4, geom="histogram",
#               #log = "x",
#               
#               binwidth=(function(x) BINWIDTH),
#               xlab = "Variant Allele Frequency (x100 %)",
#               main = paste("Hitchhikers with beneficial mutation coming in later"
#                            # #" s =", benefit,
#                            # "\n U =", u*(1-1/(1+ratio_neutral_beneficial)),  "          lifespan =",lifespan,"          population size =",
#                            # Population_size,"\n Number of Participants =",
#                            #Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#               )) + scale_y_continuous(trans = log_trans(),
#                                                                                                                         breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                                                         labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity",
#                                                                                                                                                                                               breaks = c(-24,-20,-16,-12),
#                                                                                                                                                                                               labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                                               limits =  c(lower_xlim,1)
#                                                                                                                         )
# plot7<- qplot(logq5, geom="histogram",
#               #log = "x",
#               
#               binwidth=(function(x) BINWIDTH),
#               xlab = "Variant Allele Frequency (x100 %)",
#               main = paste("Hitchhikers"
#                            # #" s =", benefit,
#                            # "\n U =", u*(1-1/(1+ratio_neutral_beneficial)),  "          lifespan =",lifespan,"          population size =",
#                            # Population_size,"\n Number of Participants =",
#                            #Number_of_Participants,"          code ran for ",code_ran_for, "mins"
#               )) + scale_y_continuous(trans = log_trans(),
#                                                                                                                         breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                                                         labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity",
#                                                                                                                                                                                               breaks = c(-9,-7,-5,-3,-1,1),
#                                                                                                                                                                                               labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                                               limits =  c(lower_xlim,1)
#                                                                                                                         )+stat_function( 
#                                                                                                                           fun = function(x){
#                                                                                                                             #1/x
#                                                                                                                             #Population_size*u*exp(-exp(x)/n_beneficial)
#                                                                                                                             BINWIDTH*Number_of_Participants*u_neu*n_hitchhiker*exp(-exp(x+log(Population_size))/n_hitchhiker)*(1/(exp(x+log(Population_size))*benefit)+benefit) 
#                                                                                                                                   #x is transformed to (x+log(Population_size) because of the normalization in the x axis
#                                                                                                                                   #no need to multiply by the population size of cells
#                                                                                                                           },
#                                                                                                                           xlim =  c(lower_xlim,-6.5)
#                                                                                                                           #, 
#                                                                                                                           #args = c(mean = mean, sd = sd, n = n, bw = binwidth)
#                                                                                                                         )
# 
# 
# plot8<- qplot(logtotal, geom="histogram",
#       #log = "x",
#       
#       binwidth=(function(x) BINWIDTH),
#       xlab = "Variant Allele Frequency (x100 %)",
#       main = paste("All mutations"
#         # #" s =", benefit,
#         # "\n U =", u,  "           lifespan =",lifespan,"          population size =",
#         # Population_size,"\n Number of Participants =",
#         # Number_of_Participants,"          code ran for ",code_ran_for, "hours"
#         )) + scale_y_continuous(trans = log_trans(), 
#                                                                                               breaks = trans_breaks("log", function(x) exp(x)),
#                                                                                               labels = trans_format("log", math_format(e^.x))) + scale_x_continuous(trans = "identity", 
#                                                                                                                                                                     breaks = c(-9,-7,-5,-3,-1,1),
#                                                                                                                                                                     labels =  trans_format("identity", math_format(e^.x)),
#                                                                                                                                                                     limits =  c(lower_xlim,1)
#                                                                                               )
# plot9<- qplot(
#               #log = "x",
# 
# 
#               main = paste(
# 
#                            "\n U =", u_ben,  "           s =", benefit,"          lifespan =",lifespan,"          population size =",
#                            Population_size,"\n Number of Participants =",
#                            Number_of_Participants,"          code ran for ",code_ran_for,"hours"))
# 
# 
# #For the curves overlaid, function(x): x is log f=l (manually transformed when being input), 
# #x is transformed to x+log(Population_size) because of the normalization over population size
# #y is transformed to y*BINWIDTH because of binsize normalization NB: given y=fun=f(x), the code plots log(y) 
# grid.arrange(plot1,plot3,plot7,nrow = 3)
# 
# grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,nrow = 5)
# 
# beep()
# 
# ######################################################################################################################
# 
# 
# #save the data 
# 
# 
# # m1 <- matrix(logp, ncol=1, byrow=TRUE)
# # d1 <- as.data.frame(m1, stringsAsFactors=FALSE)                 
# # #ggplot(data=d1,aes(d1$v1))+ geom_histogram(stat = "bin", binwidth = NULL, bins = NULL) +geom_bar(stat = "identity")
# # # (m<-
# # hist(logp,
# #     main="Distribution of mutation frequencies",
# #    xlab="mutation frequencies",
# #  border="blue",
# #   col="blue",
# #  xlim=c(0,10),
# #  #ylim=c(0,200),
# #   breaks=10000
# #  )# number of breaks from 0 to N
# # # #)
# # 
# # 
# # eq = function(x){exp(6)} 
# # lines( c(lower_xlim:1),eq(lower_xlim:1), type='l')
# # eq = function(x){exp(-s*feedingpopulation*U*x )/x}   #the only fitting comes from scaling in the y axis
# # lines(c(-500:700),eq(-500:700), type='l') #compare with theory
# #legend(300, 700, legend=c("Simulated data", "Theoretical prediction"),
# #       col=c("green", "black"), lty=1:1, cex=0.5)

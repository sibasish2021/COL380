Number of Approaches=4

1:Idea: In my first attempt i tried to write just a sequential code which would take the blocks and contruct the entire matrix and multiply to find outputs in the resulting block.
1:Why attempted:I attempted this so as to come up with a base could which always gives the correct output so that i can build up on it.
1:Results(Specify speedup): It took about 4mins to run this code on input2
1:Drawback:This used a lot of space and also took unnnecessary time due to not using cache coherence sufficiently and use of compplex data structures.

2:Idea: Here I tried to implement a sequential code with a slighltly better algorithm using blocks to find the answer of other blocks and setting the blocks only if required.
2:Why attempted:I attempted this to write a code which uses an optimised algorithm which uses the most of the fact that the matrix is sparse
2:Results(Specify speedup):This reduced my time from 4mins to about 2 min
2:Drawback:This does not use any parallelism.

3:Idea:Here I implemented parallelism by using openmp tasks to do the innner block multiplications as they can be divided into independent tasks.
3:Why attempted:I attempted this to indroduce parallelism in my code so as to use the fact that posrtions of my code can be run independently.
3:Results(Specify speedup):The speedup was quite great as I reduced from 2 mins to about 40 s.
3:Drawback:This uses parallelism but due to inherent delays because of use of 2d vectors some delay was introduced.

4:Idea:I modified the drawbacks of previous approach and used 1D vector in place of 2D vector as 1D vector has a better cache hit rate.
4:Why attempted:I attempted this to improve cache coherncy and exploit Spatial localty so as to improve cache hits.
4:Results(Specify speedup):The time decreased from 40s to around 25s.
4:Drawback:This still used many unnecessary computation as it iterated through the entire list of blocks.

5:Idea:Here I have utilised maps to iterate the list of blocks instead of direct iteration . This reduces the time complexity of the outer loops significantly. Also I used parallelism on several other loops to divide work among the threads.
5:Why attempted:I attempted this to reduce my loop iteration time which due to complete itaration of the vector of blocks was taking more time.
5:Results(Specify speedup):The time reduced from 25s to 20s,
5:Drawback:The drawback is that as we use dynamic structures like vector time due to memory latency might be more . This can be reduced by using arrays allocated space using malloc.

Scalability Analysis(time in ms)
Non-0 input blocks   Non-0 output blocks  2 cores       4 cores      8 cores     16 cores
2^10                       2^12           41562.35      20583.28     14075.72    11881.76                               
2^15                       2^17           61977.56      34577.72     20344.12    13166.48                                
2^20                       2^22           206510.88     114245.56    65146.52    45866.76                                    
2^25                       2^27           301339.24     165470.68    97678.24    61676.56                            


Thus as the number of cores is increasing the time of computation is decreasing following an almost linear trend. Thus this shows that my code is scalable.



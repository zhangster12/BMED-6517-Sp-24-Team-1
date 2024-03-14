Methods Comparison Library README 

Motivation:  

Central repository for papers used for technical benchmarking of algorithms in noninvasive cardiac signal processing. This can include ECG, SCG (accelerometer and gyro), PPG, BCG, PCG, etc. The goal is that we can have one copy of the code and then use it whenever we need a comparison for our current processing methods. These will get updated over time. As new relevant papers are found, you can add them to the list. Anyone with an itch to code can help out in implementing the algorithms. These can be great undergrad and new grad student projects.  

BIG IDEA: at a minimum, we should be at their level when it comes to processing these cardiac signals. If we don’t have a custom method that is better, let’s use the best methods other people have created. Let’s follow best practices and have an established signal processing pipeline in the lab GitHub repo.  

 

Adding a paper to the library (in MS Teams):  

Please follow the pattern in the Excel document, including the naming convention Author1_year. If adding the paper would create a duplicate, add another number to the end: Author1_year_2.  

 

Coding requirements:  

All the comments. The code should be heavily commented, so that anyone can follow what is happening and clearly see the links between the referenced paper and the code (as much as possible). If there are methods that are not able to be replicated, or multi-step processes where we leave out a few steps, make sure that is well documented.  

Each method should be a self-contained function with well-defined inputs and outputs.  

IF THE METHOD YOU ARE LOOKING AT IS VERY SIMILAR TO ANOTHER METHOD THAT HAS ALREADY BEEN CODED ask yourself if the information gain is worth the time it will take to code it up. No need to duplicate efforts. We’re looking for the best performance benchmarks, and ideas we can work off of in the future.  

Known limitations should be added to the code, whether they are from the paper or from empirical testing.  

Please add options to output figures to show that the function works as shown in the paper.  

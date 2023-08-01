<h1>Motor learning of the optokinetic response

 ### [Project Report](https://drive.google.com/file/d/1GXWopoBF_Tq_2QgXBT8rVVM08dhgQvpK/view?usp=sharing)
 
<h2>Description</h2>
The optokinetic response (OKR) is a series of rapid eye movements that occurs when the environment around us changes very quickly. Studying OKR allows us to get a better understanding of how the eyes work in combination with the brain to process visual information. This project is based on the paper "Modeling Memory Consolidation During Posttraining Periods in Cerebellovestibular Learning" [1], where authors developed a simple model of the cerebellar network. 

<p align="center">
  <img src="https://i.imgur.com/40qQot8.png" height="40%" width="40%" alt="Poster"/>
<p align="center">
  Cerebellar circuit visualization (Yamazaki et al., 2015). Variables defined in the report.

We reimplemented the model in MATLAB, and used it to further explore questions related to memory formation during and after OKR training. We were able to arrive at several new insights about the role of the cerebellar circuit in improving the optokinetic reflex. For example, we found that in PF-LTD deficient mice, which are unable to form long-term memories, having more short periods of training was better at increasing OKR gain than having less but longer periods of training, as seen in the figure below.


<p align="center">
  <img src="https://i.imgur.com/oquB58W.png" height="70%" width="70%" alt="Poster"/>
<p align="center">
  The effect of different training paradigms in short-term memory formation.

This might be a little counter intuitive, because PF-LTD deficiency affects long-term memory consolidation, so OKR improvements should not add up and the results for all training  paradigms should be the same. However, PF-LTD deficiency allows for the short-term increase in OKR to last for a longer amount of time. Therefore, if training sessions happen during the same day, like in our simulation, it is still possible for OKR improvements to add on top of each other. The rest of our findings are included in the attached project report.

<h2>Acknowledgements</h2>
This was a final project for BIOL 378 at CWRU: Computational Neuroscience, taught by Dr. Peter Thomas. The project was done in conjunction with my classmate Daniel Popp.

library("DiagrammeR")

grViz("
 digraph boxes_and_circles{
 graph[rankdir = TB]
 
 #creates nodes for schematic
 node [shape = circle, fixedsize = T, width = 1.75, penwidth= 1.5, color= cadetblue, peripheries=2, nodesep = 0.5,rankdir = TB]
     Process_Model[label= 'Process Model',fontsize= 16]
      Next_State[label= 'Next State', fontsize= 18]
     Updated_State[label= 'Updated State', fontsize= 16]
     
  node [shape = circle, width = 1.75, penwidth= 1.5, color= lightblue, nodesep = 0.5, peripheries=2]
      Observation_Model[label= 'Observation\n Model', fontsize= 16]
      Data[fontsize= 25]
     
 node [shape = diamond, height= 1, width= 1.5, color= firebrick,  peripheries=1]
     Likelihood [fontsize= 16]
     State_Space_Filter[label= 'State-Space\n Filter', fontsize= 16]
      
 node [shape = box, fixedsize = T, width = 1.25, penwidth= 1.5, color= yellowgreen, nodesep = 0.2]
     Life_history[label= 'Life History\n parameters']
     Process_Error[label= 'Process Error']
     Observation_Error[label= 'Observation Error']
  
 node [shape = box, height= 1.5,width = 2, penwidth= 0.8, color= indigo, nodesep = 0.2,peripheries=3]
    Non_detection[label= 'Non-detection\n parameter',fontsize= 20] 
     
 #add arrows
  subgraph{rank=same; Process_Model-> Next_State ->State_Space_Filter->Updated_State->Process_Model}

 Data ->Likelihood;
 
  subgraph{rank=same; Observation_Model -> Likelihood[arrowtail = T] 
  }
  
   subgraph{
   Observation_Model->Data[dir=back, minlen= 0.2]   }
  
   Next_State->Observation_Model
   
  subgraph{
   State_Space_Filter ->Likelihood[dir = both, minlen= 0.5] }
   
  
   Life_history ->Process_Model ;   Process_Error->Process_Model
   
 subgraph{
   rank = same;Observation_Error-> Observation_Model; Non_detection->Observation_Model
}

#adjust size of graph using 
graph [nodesep = 0.8]
#[label = 'Parameter estimate / Prior distribution', fontsize= 16]
 }     
      ") 


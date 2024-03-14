# Smart (Classes)

Supporting classes for the Smart package.

## Summary

This folder contains various enumerations which may be used with functions in the Smart package, as well as classes for creating useful objects.

## Table of Contents

1. [Classes](#classes)
2. [Enumerations](#enumerations)

## Classes

### Factor

**Path:** `Smart/Classes/@Factor`

The `Factor` class contains properties and methods for factor nodes -- characterized by their associated probability tables -- in graphical models. The `Factor` class is typically called by the [`Node`](#node) class when creating factor nodes.

| Function | Purpose |
| --- | --- |
| `Factor()` | Class constructor |
| `Factor.getSubtable()` | Get a subtable from a factor table |
| `Factor.makeDistribution()` | Make factor table into probability distribution |
| `Factor.makeFactor()` | Create a factor table |
| `Factor.multiply()` | Perform factor multiplication for two `Factor` or [`Message`](#message) objects |

### Node

**Path:** `Smart/Classes/@Node`

The `Node` class contains properties and methods for variable and factor nodes in graphical models. Both variable and factor nodes may be created using the class constructor. At the time of instantiation, the parents of the instantiated node may be provided.

```matlab
parent1 = Node("Parent_1", 'values', {Class.true, Class.false}); parent1 = parent1.define();
parent2 = Node("Parent_2", 'values', {Class.true, Class.false}); parent2 = parent2.define();
child = Node("Child", 'values', {0, 1, 2}); child = child.define();
```

| Function | Purpose |
| --- | --- |
| `Node()` | Class constructor |
| `Node.define()` | Define the node object |
| `Node.evaluate()` | Evaluate a factor or variable node given an array of [`Message`](#message) objects |
| `Node.quantize()` | Set the node value to a quantized version of the input |
| `Node.query()` | Query the value of the node |
| `Node.setConditionals()` | Set the conditional proababilities of the node given a data table |

### Message

**Path:** `Smart/Classes/@Message`

The `Message` class contains properties and methods for messages, for use in message-passing algorithms. `Message` objects are typically created by [`Graph`](#graph) when solving graphs.

| Function | Purpose |
| --- | --- |
| `Message()` | Class constructor |
| `Message.createMessage()` | Create a message using factor tables |

### Graph

**Path:** `Smart/Classes/@Graph`

The `Graph` class contains properties and methods for graph objects, which are collections of interconnected `Node` objects. This class is used to perform high-level operations such as querying and solving graphs.

| Function | Purpose |
| --- | --- |
| `Graph()` | Class constructor |
| `Graph.eliminate()` | Perform variable elimination on the graph |
| `Graph.query()` | Query nodes in the graph |
| `Graph.solve()` | Recursively solve a factor graph |

### Markov

**Path:** `Smart/Classes/@Markov`

The `Markov` class contains properties and methods for hidden Markov models (HMMs). After creating a model with the constructor, synthetic data may be generated, and hidden state inference can be performed using a naïve method or the Viterbi algorithm.

| Function | Purpose |
| --- | --- |
| `Markov()` | Class constructor |
| `Markov.generate()` | Generate synthetic data using a trained HMM |
| `Markov.infer()` | Infer hidden states using a naïve method |
| `Markov.viterbi()` | Infer the optimal sequence of hidden states which explains the data |

### TemplateSet

**Path:** `Smart/Classes/@TemplateSet`

The `TemplateSet` class creates sets of template signals, as described in the 2019 paper "A Unified Framework for Quality Indexing and Classification of Seismocardiogram Signals" in [IEEE JBHI](https://ieeexplore.ieee.org/document/8777167). As an example, use-case, consider that sets of signals are available as column matrices for several subjects. To create a template set from this data, we first format the data as a cell vector, where each cell contains a column matrix with the data from each subject. The distance metric for the SQI can be specified with the [`Index`](#index) enumeration.

```matlab
templateSet = TemplateSet({data}, <Index>, <lambda>);
```

Note that `<lambda>` is typically assigned 25 for cardiac bio-signals, and that `data` is nested in extra brackets. This is because each cell in the external bracket corresponds to data from a different class. If there is only one class of data, the external cell array thus contains one element. To create a template set from the data, the command `templateSet = templateSet.create();` may be used. To generate scores on an array of signals from another subjects, the command `scores = templateSet.score('data', {new_data});` may be used. Scores are returned in a cell array corresponding to each template. To obtain a mean score for each held-out signal, we call the function `mean_scores = templateSet.meanScore(scores)`.

For classification tasks, the dataset may be split by class, as aforementioned.

```matlab
templateSet = TemplateSet({data_c1, data_c2 ...}, <Index>, <lambda>)
```

The template set must then be characterized to perform classification tasks. To do so, call the function `templateSet = templateSet.characterize();` and `templateSet = templateSet.predict('data', {new_data_c1, new_data_c2, ...});` to obtain predictions for each cell in the vector. For additional functionality of the `TemplateSet` class, type `help TemplateSet.<function>` in the console.

| Function | Purpose |
| --- | --- |
| `TemplateSet()` | Class constructor |
| `TemplateSet.characterize()` | Train the `TemplateSet` object for subsequent prediction |
| `TemplateSet.confidence()` | Obatin the probability of prediction error versus the number of observations after prediction |
| `TemplateSet.confusion()` | Create a confusion matrix from the template set |
| `TemplateSet.create()` | Create a set of templates from the data |
| `TemplateSet.predict()` | Predict class labels with a trained `TemlpateSet` object |
| `TemplateSet.score()` | Generate SQI scores with a trained `TemplateSet` object |

### Manifold

**Path:** `Smart/Classes/@Manifold`

The `Manifold` class characterizes manifolds and extracts latent variables using the ISOMAP method. To use this class, we first initialize the `Manifold` object using `manifold = Manifold("Description")` and load the data using `manifold = manifold.createGraph(data)`. The resulting graph may be visualized with the function `manifold.plotGraph()`. To implement the ISOMAP algorithm, we must first compute the geodesic distances between all nodes in the graph using the function `manifold = manifold.shortestPath()`, followed by latent variable extraction via `manifold = manifold.scale()`, which performs multidimensional scaling.

| Function | Purpose |
| --- | --- |
| `createGraph()` | Create adjacency graph for performing ISOMAP |
| `latent()` | Function returning latent variables and corresponding original datapoints |
| `map()` | Map new datapoints to learned manifold |
| `plotGraph()` | Plot datapoints on adjacency graph |
| `scale()` | Perform multidimensional scaling using dissimilarity matrix |
| `shortestPath()` | Computing geodesic distances between points in adjacency graph |

## Enumerations

### Index

**Path:** `Cardio/Classes/Class.m`

Enumeration describing possible boolean values for `Node` objects.

| Value | Meaning |
| --- | --- |
| true | Boolean true |
| false | Boolean false |
| nil | Boolean nil |

### MessageEval

**Path:** `Cardio/Classes/MessageEval.m`

Enumeration describing methods of message evaluation.

| Value | Meaning |
| --- | --- |
| sumProduct | Sum-product algorithm |
| maxProduct | Max-product algorithm |
| none | N/A |

### MessageType

**Path:** `Cardio/Classes/MessageType.m`

Enumeration describing types of messages.

| Value | Meaning |
| --- | --- |
| varToVar | Variable node to variable node |
| varToFac | Variable node to factor node |
| facToVar | Factor node to variable node |

### NodeType

**Path:** `Cardio/Classes/NodeType.m`

Enumeration describing types of `Node` objects.

| Value | Meaning |
| --- | --- |
| Variable | Variable node |
| Factor | Factor node |

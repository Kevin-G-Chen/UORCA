In this document, I will be making notes about the WA Government's AI policy, and I will consider how my work/work needs to adhere to these policies.

# Human, social and environmental wellbeing

>Throughout their lifecycle, AI systems must benefit individuals, society, and the environment.

The primary purpose of UORCA is to assist with the interpretation of biological data, with the downstream impact of improving health outcomes for patients across WA (and nationally/internationally). I see my project meeting this guideline.

# Human-centred values

>Throughout their lifecycle, AI systems must respect human rights, diversity, and the autonomy of individuals

Hm, this one is tricky. The fundamental considerations are that "how could my AI system affect human rights?" ... and I feel thinking about this more, I don't think there is any realistic universe where this happens.

# Fairness

>Are the models trained and tested on relevant, accurate, and generalisable datasets and is the AI system deployed by users trained to implement them responsibly to manage and mitigate bias

1. There is explicit mention of ensuring I am testing (in my case, benchmarking) on good datasets. This in turn also seems to imply (not that I wouldn't do so) that I should make a record of the datasets that I am performing the benchmarking on
2. A bit tricky to decipher what is meant by the second part. It seems to be about mitigating biases - therefore, what this means is that I should take into account where might biases appear in my results, and what steps should I take to ensure minimise bias.

# Privacy protection and security

>Compliance with appropriate data policies and legislation, for example, the forthcoming WA Privacy & Responsible Information Sharing (PRIS) and the Commonwealth Privacy Act 1988

I will need to consider privacy/data security measures. Will the workflow be resilient to attacks?

Immediately I'm thinking of the fact that, since I am downloading files, I need to consider if this will compromise any security.

# Reliability and safety

>Throughout their lifecycle, AI systems myst reliably operate in accordance with their intended purpose.

One of the points that is raised here is that AI systems must be monitored - this wasn't something that I had considered before.

There is a point about safety risks - I don't expect that my system would pose safety risks, but I should think about if there is a formal way to determine the potential safety risks at each step.

They make a point about making clear the connection between data and inferences drawn from that data. This must be "sound" and "assessed in an ongoing manner"
... but I'm not sure how I will manage that.

# Transparency, explainability, and contestability

>So affected stakeholders can know how the AI model reached its decision

People must have access to an efficient and transparent review mechanism - so what constitutes a "review mechanism"? It seems this refers to knowing when AI is involved.

There is also a point made about having sufficient access/understanding of the algorithm. I internally interpret this as meaning I should ensure my code is openly available?

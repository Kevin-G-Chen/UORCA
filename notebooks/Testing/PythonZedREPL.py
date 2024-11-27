# %% Cell 1
# Just testing the Python REPL in Zed.

from openai import OpenAI
import os
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import numpy as np
# %% Cell 2
load_dotenv('~/Documents/UORCA/.env')
openai_api_key = os.getenv('OPENAI_API_KEY')
# %% Cell 3
# Test the OpenAI API key is working

client = OpenAI(
  api_key=openai_api_key,  # this is also the default, it can be omitted
)
chat_completion = client.chat.completions.create(
    messages=[
        {
            "role": "user",
            "content": "What is the capital city of Victoria?",
        }
    ],
    model="gpt-4o-mini"
)

result = chat_completion.choices[0].message.content
print(result)  # Testing
# Basic plot using matplotlib

x = np.linspace(0, 10, 100)
y = np.sin(x)

plt.plot(x, y)

plt.title('Simple Sine Wave')
plt.xlabel('x')
plt.ylabel('sin(x)')
plt.grid(True)
plt.show()

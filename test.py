print('Hello Word !')
import pandas as pd
data = {
    'Name': ['Alice', 'Bob', 'Charlie'],
    'Age': [25, 30, 35],
    'City': ['New York', 'San Francisco', 'Los Angeles']
}

# Create a DataFrame
df = pd.DataFrame(data)

# Display the DataFrame
print(df)
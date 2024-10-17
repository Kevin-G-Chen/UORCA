import matplotlib.pyplot as plt

def fibonacci(n):
    """Generate Fibonacci sequence up to n terms."""
    sequence = [0, 1]
    if n <= 0:
        return []
    elif n == 1:
        return [0]
    elif n == 2:
        return sequence

    for i in range(2, n):
        next_term = sequence[-1] + sequence[-2]
        sequence.append(next_term)

    return sequence

def get_user_input():
    """Get number of Fibonacci terms from the user."""
    while True:
        try:
            n = int(input("Enter the number of Fibonacci terms to generate: "))
            if n <= 0:
                print("Please enter a positive integer.")
            else:
                return n
        except ValueError:
            print("Invalid input. Please enter a valid integer.")

def plot_fibonacci(sequence):
    """Plot the Fibonacci sequence."""
    plt.figure(figsize=(10, 6))
    plt.plot(sequence, marker="o", linestyle="-", color="b", label="Fibonacci Sequence")
    plt.title(f"Fibonacci Sequence (First {len(sequence)} Terms)")
    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.grid(True)
    plt.legend()
    plt.show()

def main():
    print("Fibonacci Sequence Generator")
    
    # Get the number of terms from the user
    n = get_user_input()
    
    # Generate the Fibonacci sequence
    fib_sequence = fibonacci(n)
    
    # Print the sequence
    print(f"Fibonacci Sequence (First {n} Terms):")
    print(fib_sequence)
    
    # Plot the sequence
    plot_fibonacci(fib_sequence)

if __name__ == "__main__":
    main()

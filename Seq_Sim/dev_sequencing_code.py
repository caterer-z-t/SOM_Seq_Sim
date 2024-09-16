from utils import *
import sys


def main():
    print(optimal_subplot_layout(9))
    print(optimal_subplot_layout(10))
    print(optimal_subplot_layout(11))
    print(optimal_subplot_layout(12))
    input_file = sys.argv[1]

    df = pd.read_csv(input_file, sep="\t")

    plot_numerical_distributions(metadata=df)

    plot_categorical_distributions(metadata=df)


if __name__ == "__main__":
    main()

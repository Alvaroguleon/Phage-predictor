import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import pickle
import argparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, f1_score,  log_loss, mean_squared_error, balanced_accuracy_score
from math import sqrt

def matthews(y_true, y_pred):
    """
    Calculate the Matthews Correlation Coefficient and other metrics.
    """
    if type(y_true) == pd.Series:
        y_true = y_true.values

    P = len([x for x in y_true if x == 1])
    N = len([x for x in y_true if x == 0])

    Tp, Fp = 0, 0
    for i in range(len(y_true)):
        if y_true[i] == 1 and y_pred[i] == 1: Tp += 1
        elif y_true[i] == 0 and y_pred[i] == 1: Fp += 1

    Tn = N - Fp
    Fn = P - Tp

    try:
        mcc = (Tp * Tn - Fp * Fn) / sqrt(
            (Tn + Fn) * (Tn + Fp) * (Tp + Fn) * (Tp + Fp))
    except ZeroDivisionError:
        mcc = 0

    return (mcc, f" \n \
    P: {P:_} \n \
    Tp: {Tp:_} \n \
    Fp: {Fp:_} \n \
    N: {N:_} \n \
    Tn: {Tn:_} \n \
    Fn: {Fn:_}")


#TODO: change output name of reports and models to reflect current day, instead of being fixed names
parser = argparse.ArgumentParser(description="Train the model using processed data.")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
# optional.add_argument("-d", "--data", type=str, required=True,
#                     help="Path to the data to use in the model (csv).")
optional.add_argument("-d","--data", dest="data", action='store', default='data/processed/model_data.csv', help="Path to the data to use in the model (default: data/processed/model_csv).")
optional.add_argument("-o","--output_directory", dest="output", action='store', default='models/random_forest.pkl', help="Output path for the model (default: models/random_forest.pkl).")
optional.add_argument("-r","--report", dest="report", action='store', default='reports/training_random_forest.txt', help="Output path for the model (default: reports/training_random_forest.txt).")

args = parser.parse_args()

print("--- Initiating training of the model ---")

df = pd.read_csv(args.data, index_col=0)

# Features (independent variables)
features = ['genome_length', 'jumbophage', 'gc_%',
       'trna_count', 'cds_number', 'coding_capacity', 'positive_strand_%',
       'negative_strand_%', 'molecule_type_ss-DNA', 'molecule_type_DNA',
       'molecule_type_RNA', 'molecule_type_ss-RNA', 'topology_circular','topology_linear']

# Target variable (dependent variable)
target = 'staining'

# Extract features and target
X = df[features]
y = df[target]

#TODO: Change it to introduce homology partitioning

# Split the data into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a logistic regression model
rf = RandomForestClassifier()
rf.fit(X_train, y_train)

print("Saving model...")
# save model
pickle.dump(rf, open(args.output, "wb"))



# Predict probabilities on the test set
y_pred_proba = rf.predict_proba(X_test)[:, 1]  # Probabilities for the positive class

# Predict class labels for F1 score calculation
y_pred = rf.predict(X_test)

# Encode labels for mean squared error calculation
label_encoder = LabelEncoder()
y_test_encoded = label_encoder.fit_transform(y_test)
y_pred_encoded = label_encoder.transform(y_pred)


# Writing report --------
# Calculate F1 score
positive_label = 'positive'
f1score = f1_score(y_test, y_pred, pos_label=positive_label)

# Calculate log loss using predicted probabilities
logloss = log_loss(y_test, y_pred_proba)

# Calculate mean squared error using encoded labels
mse = mean_squared_error(y_test_encoded, y_pred_encoded)

# Calculate accuracy
accuracy = accuracy_score(y_test_encoded, y_pred_encoded)

# Calculate balanced accuracy
balanced_accuracy = balanced_accuracy_score(y_test_encoded, y_pred_encoded)

# Calculate Matthews Correlation Coefficient (MCC)
mcc, mcc_details = matthews(y_test_encoded, y_pred_encoded)

# Generate the classification report
classification_report_str = classification_report(y_test_encoded, y_pred_encoded)


report_str = f"""
Classification Report:
{classification_report_str}

Accuracy: {accuracy}
Balanced Accuracy: {balanced_accuracy}
Log Loss: {logloss}
Mean Squared Error: {mse}
F1 Score: {f1score}
Matthews Correlation Coefficient: {mcc}
Matthews Correlation Coefficient Details: {mcc_details}
"""

# Save the report to a text file
with open(args.report, "w") as report_file:
    report_file.write(report_str)

print(report_str) 




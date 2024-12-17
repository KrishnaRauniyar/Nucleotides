import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, StandardScaler
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
from tensorflow.keras.regularizers import l2
import argparse

def drugClassifier(file_path):
    def parse_protein_name(protein_name):
        parts = protein_name.split('_')
        return parts[0], parts[1], parts[2], parts[3]
    
    df = pd.read_csv(file_path)

    # Parsing protein names and adding them as separate columns
    df[['protein_id', 'chain_id', 'residue_number', 'residue_type']] = df['Protein Name'].apply(lambda x: pd.Series(parse_protein_name(x)))

    # Dropping the original 'Protein Name' column
    df.drop(columns=['Protein Name'], inplace=True)

    # Encoding categorical variables
    label_encoders = {}
    for column in ['protein_id', 'chain_id', 'residue_type']:
        le = LabelEncoder()
        df[column] = le.fit_transform(df[column])
        label_encoders[column] = le

    # Printing the total number of occurrences of each residue_type with their labels
    residue_type_counts = df['residue_type'].value_counts()
    residue_type_labels = label_encoders['residue_type'].inverse_transform(residue_type_counts.index)
    residue_type_count_dict = dict(zip(residue_type_labels, residue_type_counts))
    print("Total number of occurrences of each residue_type:")
    for label, count in residue_type_count_dict.items():
        print(f"{label}: {count}")

    # Defining features and target
    X = df.drop(columns=['residue_number', 'residue_type'])
    y = df['residue_type']

    # Standardizing the feature matrix
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Splitting the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    # Defining and compile the neural network model
    model = Sequential([
        Dense(128, input_dim=X_train.shape[1], activation='relu'),
        Dense(64, activation='relu', kernel_regularizer=l2(0.01)),
        Dense(32, activation='relu', kernel_regularizer=l2(0.01)),
        Dense(16, activation='relu', kernel_regularizer=l2(0.01)),
        Dense(8, activation='relu', kernel_regularizer=l2(0.01)),
        Dense(len(label_encoders['residue_type'].classes_), activation='softmax')
    ])
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])

    # Adding early stopping
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

    # Training the model
    history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_data=(X_test, y_test), callbacks=[early_stopping])

    # Evaluating the model
    loss, accuracy = model.evaluate(X_test, y_test)
    print("-----------------------------------------")
    print(f'Accuracy: {accuracy}')

    # Plotting training and validation loss
    history_dict = history.history
    loss = history_dict['loss']
    val_loss = history_dict['val_loss']
    plt.plot(loss, label='Training Loss')
    plt.plot(val_loss, label='Validation Loss')
    plt.title('Training and Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    # plt.ylim([0, max(max(loss), max(val_loss))])
    plt.savefig("loss_plot.png")
    plt.close()

    # Plotting training and validation accuracy
    accuracy = history_dict['accuracy']
    val_accuracy = history_dict['val_accuracy']
    plt.plot(accuracy, label='Training Accuracy')
    plt.plot(val_accuracy, label='Validation Accuracy')
    plt.title('Training and Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend(loc = 'lower left')
    # plt.ylim([0.7, 1])
    plt.savefig("accuracy_plot.png")
    plt.close()

    # Making predictions and generating confusion matrix and classification report
    y_pred = model.predict(X_test)
    y_pred_classes = y_pred.argmax(axis=1)
    class_labels = label_encoders['residue_type'].classes_
    cm = confusion_matrix(y_test, y_pred_classes, labels=np.arange(len(class_labels)))
    cr = classification_report(y_test, y_pred_classes, labels=np.arange(len(class_labels)), target_names=class_labels)
    print(cm)
    print(cr)

    # Creating confusion matrix display
    cmd = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=class_labels)
    # Plotting the confusion matrix
    fig, ax = plt.subplots(figsize=(25, 25))  # Adjusting the figsize to better fit the labels
    cmd.plot(cmap=plt.cm.Blues, ax=ax)
    plt.title('Confusion Matrix')
    plt.xticks(np.arange(len(class_labels)), class_labels, rotation=90, fontsize=8)
    plt.yticks(np.arange(len(class_labels)), class_labels, fontsize=8)
    plt.savefig("confusion_matrix.png")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser("DNN Model")
    parser.add_argument('-p','--input_path', type=str, required=True, help="CSV input path")
    args = parser.parse_args()
    drugClassifier(args.input_path)

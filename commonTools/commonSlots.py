# simple and general slots called by gui.py

import os
from PySide2.QtWidgets import QLineEdit, QFileDialog, QMessageBox

def openFileDialog(fileType, lineEdit):
    ''' return a openFile function that opens special type of files
        the openFile dialog is connected with a lineEdit widget
        inputs: 
            fileType (string): the file type for opening
            lineEdit (QLineEdit): the QLineEdit that connects to the QFileDialog '''
    def openFile():
        fileName, _ = QFileDialog.getOpenFileName(
            None,
            f'Choose {fileType} file',
            os.getcwd(),
            f'All Files (*)'
        )
        lineEdit.setText(fileName)
    return openFile

def openFilesDialog(fileType, listWidget):
    ''' return a openFile function that opens a series of files
        the openFile dialog is connected with a QListWidget
        inputs: 
            fileType (string): the file type for opening
            listWidget (QListWidget): the QListWidget that connects to the QFileDialog '''
    def openFiles():

        fileNames, _ = QFileDialog.getOpenFileNames(
            None,
            f'Choose {fileType} files',
            os.getcwd(),
            f'All Files (*)'
        )
        listWidget.addItems(fileNames)
    return openFiles

# simple and general slots called by gui.py

from PySide2.QtWidgets import QLineEdit, QFileDialog, QMessageBox

def openFileDialog(fileType, lineEdit):
    """return a openFile function that opens special type of files
       the openFile dialog is connected with a lineEdit widget

    Args:
        fileType (str): the file type for opening
        lineEdit (QLineEdit): the QLineEdit that connects to the QFileDialog

    Returns:
        function obj: slot function opening special type of files
    """    
    
    def openFile():
        fileName, _ = QFileDialog.getOpenFileName(
            None,
            f'Choose {fileType} file',
            '',
            f'All Files (*)'
        )
        lineEdit.setText(fileName)
    return openFile

def openFilesDialog(fileType, listWidget):
    """return a openFile function that opens a series of files
       the openFile dialog is connected with a QListWidget

    Args:
        fileType (str): the file type for opening
        listWidget (QListWidget): the QListWidget that connects to the QFileDialog

    Returns:
        function obj: slot function opening a series of files
    """    

    def openFiles():

        fileNames, _ = QFileDialog.getOpenFileNames(
            None,
            f'Choose {fileType} files',
            '',
            f'All Files (*)'
        )
        listWidget.addItems(fileNames)
    return openFiles

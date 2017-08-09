import QtQuick 2.6
import QtQuick.Dialogs 1.2
import QtQuick.Controls 1.5

ApplicationWindow {
    visible: true
    width: 640
    height: 480
    title: qsTr("Demo App")

    function openFile(fileUrl) {
        var request = new XMLHttpRequest();
        request.open("GET", fileUrl, false);
        request.send(null);
        return request.responseText;
    }

    function saveFile(fileUrl, text) {
        var request = new XMLHttpRequest();
        request.open("PUT", fileUrl, false);
        request.send(text);
        return request.status;
    }

    FileDialog {
        id: openFileDialog
        nameFilters: ["Text files (*.txt)", "All files (*)"]
        onAccepted: textEdit.text = openFile(openFileDialog.fileUrl)
    }

    FileDialog {
        id: saveFileDialog
        selectExisting: false
        nameFilters: ["Text files (*.txt)", "All files (*)"]
        onAccepted: saveFile(saveFileDialog.fileUrl, textEdit.text)
    }

    menuBar: MenuBar {
        Menu {
            title: qsTr("File")
            MenuItem {
                text: qsTr("&Open")
                onTriggered: openFileDialog.open()
            }
            MenuItem {
                text: qsTr("&Save")
                onTriggered: saveFileDialog.open()
            }
            MenuItem {
                text: qsTr("Exit")
                onTriggered: Qt.quit();
            }
        }
    }

    TextArea {
        id: textEdit
        anchors.fill: parent
        text:
            "Lorem ipsum dolor sit amet, consectetur adipisicing elit, " +
            "sed do eiusmod tempor incididunt ut labore et dolore magna " +
            "aliqua. Ut enim ad minim veniam, quis nostrud exercitation " +
            "ullamco laboris nisi ut aliquip ex ea commodo cosnsequat. ";
    }
}

import QtQuick 2.0

Item {
    width: 100; height: 30;
    property string name: ""
    property string defaultIn: ""
    Text {
        id: text
        height: parent.height
        width: parent.width / 2
        text: qsTr(name)
        font.pixelSize: 15
        anchors.fill: left
    }
    Item {
        //property alias text: textInput.text
        height: parent.height
        width: parent.width/2
        anchors.fill: right
        BorderImage {

            id: borderim
            source: "images/lineedit.png"
            width: parent.width; height: parent.height
            border.left: 10; border.top: 10
            border.right: 10; border.bottom: 10

        }
        TextInput {
            id: textInput
            width: parent.width-10
            height: parent.height-5
            text: qsTr(defaultIn)
            selectionColor: "green"
            font.bold: true
            font.pixelSize: 12
            maximumLength: 16
            anchors.centerIn: parent
            focus: true
        }
    }
}


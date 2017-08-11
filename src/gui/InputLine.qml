import QtQuick 2.0

Rectangle {
    id: root
    property string name: ""
    property string defaultIn: ""
    property int fontSize: 16
    width: parent.width
    height: 36
    color: "#cdcdcd"
    border.color: "black"
    border.width: 2
    Rectangle {
        id: label
        color: parent.color
        width: parent.width/4
        height: parent.height
        Text {
            anchors.verticalCenter: parent
         //   anchors.left: parent.left
         //   anchors.bottom: parent.bottom
            anchors.margins: 10
            font.pixelSize: fontSize
            text: qsTr(name)
        }
    }
    Rectangle {
        color: "white"
        radius: 5
        border.color: "black"
        border.width: 2
        anchors.left: label.right
        width: parent.width/4*3
        height: parent.height

        TextInput {
            id: textInput
            color: "black"
            text: qsTr(defaultIn)
            //font.bold: true
            font.pixelSize: fontSize
            maximumLength: 50
            //anchors.verticalCenter: parent
            anchors.left: parent.left
            anchors.margins: 10
        }
    }

}

/*
Rectangle {
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
*/

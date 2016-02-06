#include "notarypage.h"
#include "ui_notarypage.h"
#include "uint256.h"
#include "walletmodel.h"
#include "util.h"
#include "openssl/sha.h"
#include "guiutil.h"

#include <QApplication>
#include <QClipboard>
#include <QFileDialog>
#include <QMessageBox>
#include <QMenu>

NotaryPage::NotaryPage(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::NotaryPage),
    model(0)
{
    ui->setupUi(this);
    ui->tableWidget->setColumnCount(2);
    ui->tableWidget->setRowCount(0);

    // Column headers
    QString txHeader = tr("Transaction");
    QString blockHeader = tr("Block");
    QStringList tableHeaders = (QStringList() << txHeader << blockHeader);
    ui->tableWidget->setHorizontalHeaderLabels(tableHeaders);

    // Column widths
#if QT_VERSION < 0x050000
    ui->tableWidget->horizontalHeader()->setResizeMode(0, QHeaderView::Stretch);
#else
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
#endif
    ui->tableWidget->horizontalHeader()->setStretchLastSection(false);

    // Context menu
    QAction *copyTxAction = new QAction(tr("Copy To Clipboard"), this);
    contextMenu = new QMenu();
    contextMenu->addAction(copyTxAction);
    ui->tableWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(copyTxAction, SIGNAL(triggered()), this, SLOT(onCopyTxID()));
    connect(ui->tableWidget, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextualMenu(QPoint)));

}

NotaryPage::~NotaryPage()
{
    delete ui;
}

void NotaryPage::setSearchResults(std::vector<std::pair<std::string, int> > txResults)
{
    ui->searchNotaryButton->setEnabled(true);
    ui->tableWidget->setRowCount(txResults.size());

    int rowNumber = 0;
    for (std::vector<std::pair<std::string, int> >::iterator iter = txResults.begin(); iter != txResults.end(); ++iter)
    {
        QString hashQ = QString::fromStdString(iter->first);
        QString blockQ = QString::number(iter->second);

        QTableWidgetItem *txItem = new QTableWidgetItem(hashQ);
        QTableWidgetItem *blockItem = new QTableWidgetItem(blockQ);
        txItem->setFlags(txItem->flags() ^ Qt::ItemIsEditable);
        blockItem->setFlags(blockItem->flags() ^ Qt::ItemIsEditable);
        ui->tableWidget->setItem(rowNumber, 0, txItem);
        ui->tableWidget->setItem(rowNumber, 1, blockItem);

        rowNumber++;
    }
}

void NotaryPage::showNotaryTxResult(std::string txID, std::string txError)
{
    if (txError == "") {
        std::string txSentMsg = "Notary transaction sent successfully: " + txID;
        QMessageBox::information(this, tr("Send Notary Tx"),
            tr(txSentMsg.c_str()),
            QMessageBox::Ok, QMessageBox::Ok);
    } else {
        QMessageBox::warning(this, tr("Send Notary Tx"),
            tr(txError.c_str()),
            QMessageBox::Ok, QMessageBox::Ok);
    }
}

void NotaryPage::setNotaryFileName(QString fileName)
{
    ui->sendNotaryEntry->setText(fileName);
}

void NotaryPage::on_searchNotaryButton_clicked()
{
    std::string notaryID = ui->notaryIDEdit->text().toStdString();
    if (!(IsHex(notaryID) && notaryID.length() == 64)) {
        ui->notaryIDEdit->setValid(false);
        return;
    }

    uint256 hash(notaryID);
    ui->searchNotaryButton->setEnabled(false);
    model->searchNotaryTx(hash);
}

void NotaryPage::setModel(WalletModel *model)
{
    this->model = model;
    connect(this->model, SIGNAL(notarySearchComplete(std::vector<std::pair<std::string, int> >)), ui->tableWidget, SLOT(clearContents()));
    connect(this->model, SIGNAL(notarySearchComplete(std::vector<std::pair<std::string, int> >)), this, SLOT(setSearchResults(std::vector<std::pair<std::string, int> >)));
    connect(this->model, SIGNAL(notaryTxSent(std::string, std::string)), this, SLOT(showNotaryTxResult(std::string, std::string)));
}

void NotaryPage::on_selectFileButton_clicked()
{
    QString fileName;
    QFileDialog dlg(this);
    dlg.setFileMode(QFileDialog::ExistingFile);

    if (dlg.exec())
    {
        fileName = dlg.selectedFiles()[0];
        setNotaryFileName(fileName);
        ui->calcNotaryIDbutton->animateClick();
    }
}

void NotaryPage::on_sendNotaryButton_clicked()
{
    std::string fileHash = ui->notaryIDEdit->text().toStdString();
    if (!(IsHex(fileHash) && fileHash.length() == 64)) {
        ui->notaryIDEdit->setValid(false);
        return;
    }

    // Make sure wallet is unlocked
    WalletModel::UnlockContext ctx(model->requestUnlock());
    if (!ctx.isValid()) {
        return;
    }

    model->sendNotaryTx(fileHash);
}

std::string NotaryPage::hashFile(std::string fileName)
{
    unsigned char hashSha[SHA256_DIGEST_LENGTH];
    FILE* file = fopen(fileName.c_str(), "rb");
    if (file == NULL) {
        return "";
    }

    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    char buffer[4096];
    int bytesRead = 0;
    while ((bytesRead = fread(buffer, 1, 4096, file)) != 0) {
        SHA256_Update(&sha256, buffer, bytesRead);
    }
    SHA256_Final(hashSha, &sha256);
    std::string notaryID = HashToString(hashSha, SHA256_DIGEST_LENGTH);
    LogPrintf("NotaryPage::hashFile: hash of file: %s\n", notaryID);

    fclose(file);

    return notaryID;
}

void NotaryPage::contextualMenu(const QPoint &point)
{
    QModelIndex index = ui->tableWidget->indexAt(point);
    if (index.isValid())
    {
        contextMenu->exec(QCursor::pos());
    }
}

void NotaryPage::onCopyTxID()
{
    QString txID = ui->tableWidget->selectedItems().at(0)->text();
    if (txID.length() > 0) {
        QApplication::clipboard()->setText(txID);
    }
}

void NotaryPage::calculateNotaryID()
{
    std::string fileName = ui->sendNotaryEntry->text().toStdString();
    std::string fileHash = hashFile(fileName);
    // Warn if file is NULL
    if (fileHash == "") {
        QMessageBox::warning(this, tr("Notarize File"),
            tr("Unable to open file for hashing."),
            QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    ui->notaryIDEdit->setText(QString::fromStdString(fileHash));
}

void NotaryPage::on_calcNotaryIDbutton_clicked()
{
    calculateNotaryID();
}

#include "clamourpage.h"
#include "ui_clamourpage.h"
#include "openssl/sha.h"
#include "clamspeech.h"
#include "main.h"
#include "util.h"
#include "walletmodel.h"
#include "clamourpetitionmodel.h"
#include "clamoursupportmodel.h"

#include <QDebug>
#include <QMessageBox>
#include <QMenu>

ClamourPage::ClamourPage(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ClamourPage),
    model(0)
{
    ui->setupUi(this);
    ui->createPetitionButton->setEnabled(false);
    ui->setVoteCheckBox->setEnabled(false);

    petitionModel = new ClamourPetitionModel(this);
    ui->searchClamourView->setModel(petitionModel);

    // Context menu for petition support view
    QAction *searchPetitionIDAction = new QAction(tr("Search for petition"), this);
    petitionViewContextMenu = new QMenu();
    petitionViewContextMenu->addAction(searchPetitionIDAction);
    connect(searchPetitionIDAction, SIGNAL(triggered()), this, SLOT(searchHighlightedPetition()));

    supportModel = new ClamourSupportModel(this);
    ui->petitionSupportView->setModel(supportModel);
    connect(ui->petitionSupportView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(searchHighlightedPetition()));
}

ClamourPage::~ClamourPage()
{
    delete ui;
}

// Calculate notary ID when text changes.
void ClamourPage::on_createPetitionEdit_textChanged()
{
    std::string petitionText(ui->createPetitionEdit->toPlainText().toStdString());
    if (petitionText.length() == 0)
    {
        ui->petitionIDEdit->clear();
        ui->createPetitionButton->setEnabled(false);
        ui->setVoteCheckBox->setEnabled(false);
        return;
    }
    ui->createPetitionButton->setEnabled(true);
    ui->setVoteCheckBox->setEnabled(true);
    std::string petitionHash(StrToSHA256(petitionText));
    ui->petitionIDEdit->setText(QString::fromStdString(petitionHash));
}

// Create a tx that creates a petitition
void ClamourPage::on_createPetitionButton_clicked()
{
    std::string petitionHash(ui->petitionIDEdit->text().toStdString());

    WalletModel::UnlockContext ctx(model->requestUnlock());
    if (!ctx.isValid()) {
        return;
    }

    model->sendClamourTx(petitionHash);

    if (ui->setVoteCheckBox->isChecked())
    {
        strDefaultStakeSpeech = "clamour " + petitionHash.substr(0, 8);
        clamourClamSpeech.push_back(strDefaultStakeSpeech);
        qDebug() << "saving clamour petitions";
        if ( !SaveClamourClamSpeech() )
            qDebug() << "Clamour CLAMSpeech petitions FAILED to save!";
        loadVotes();
    }
}

void ClamourPage::on_setVotesButton_clicked()
{
    saveVotes();
}

void ClamourPage::loadVotes()
{
    QStringList list;
    for (std::vector<std::string>::iterator it = clamourClamSpeech.begin(); it != clamourClamSpeech.end(); ++it)
    {
        list.append(QString::fromStdString(*it).mid(8));
    }
    ui->votesEdit->setPlainText(list.join("\n"));
}

void ClamourPage::saveVotes()
{
    QStringList list = ui->votesEdit->toPlainText().replace("\n", ",").replace(" ", ",").split(',', QString::SkipEmptyParts);
    std::vector<std::string> newSpeeches;
    clamourClamSpeech.clear();

    if (list.length() > 0)
    {
        newSpeeches.push_back("clamour");
        foreach ( const QString &strLine, list )
            if ( !strLine.isEmpty() && strLine.length() >= 8 && IsHex(strLine.toStdString()) )
            {
                // Create new string if necessary
                if (newSpeeches.back().length() > MAX_TX_COMMENT_LEN - 9)
                {
                    newSpeeches.push_back("clamour");
                }
                std::string &newSpeech = newSpeeches.back();
                newSpeech = newSpeech + " " + strLine.trimmed().left(8).toStdString();
            }


        for (std::vector<std::string>::iterator it = newSpeeches.begin(); it != newSpeeches.end(); ++it)
        {
            clamourClamSpeech.push_back(*it);
        }
    }

    // save new speech
    qDebug() << "saving clamour petitions";
    if ( !SaveClamourClamSpeech() )
        qDebug() << "Clamour CLAMSpeech petitions FAILED to save!";

    loadVotes();
}

void ClamourPage::showClamourTxResult(std::string txID, std::string txError)
{
    if (txError == "") {
        std::string txSentMsg = "Clamour petition created successfully: " + txID;
        QMessageBox::information(this, tr("Create Clamour Petition"),
            tr(txSentMsg.c_str()),
            QMessageBox::Ok, QMessageBox::Ok);
        ui->createPetitionButton->setEnabled(false);
        ui->setVoteCheckBox->setEnabled(false);
    } else {
        QMessageBox::warning(this, tr("Create Clamour Petition"),
            tr(txError.c_str()),
            QMessageBox::Ok, QMessageBox::Ok);
    }
}

void ClamourPage::setClamourSearchResults(CClamour *pResult)
{
    if (!pResult)
    {
        LogPrintf("No clamour results.\n");
        QMessageBox::warning(this, tr("Clamour Search"),
            tr("No clamour petition found."),
            QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    petitionModel->setPetition(pResult);
}

bool petitionPairSort(std::pair<std::string, int> i, std::pair<std::string, int> j) { return i.second > j.second; }

void ClamourPage::showPetitionSupport(std::map<string, int> mapSupport)
{
    std::vector<std::pair<std::string, int> > support;
    for (std::map<std::string, int>::iterator it = mapSupport.begin(); it != mapSupport.end(); ++it) {
        std::pair<std::string, int> petitionSupport = std::make_pair(it->first, it->second);
        support.push_back(petitionSupport);
    }
    std::sort(support.begin(), support.end(), petitionPairSort);
    supportModel->setSupport(support);
}

void ClamourPage::setModel(WalletModel *model)
{
    this->model = model;
    connect(this->model, SIGNAL(clamourTxSent(std::string, std::string)), this, SLOT(showClamourTxResult(std::string, std::string)));
    connect(this->model, SIGNAL(clamourSearchComplete(CClamour*)), this, SLOT(setClamourSearchResults(CClamour*)));
    connect(this->model, SIGNAL(petitionSupportRetrieved(std::map<std::string,int>)), this, SLOT(showPetitionSupport(std::map<std::string,int>)));
    loadVotes();
}

void ClamourPage::on_searchClamourButton_clicked()
{
    std::string pid(ui->searchClamourEdit->text().toStdString());
    if (!(IsHex(pid) && pid.length() == 8)) {
        ui->searchClamourEdit->setValid(false);
        return;
    }
    petitionModel->clear();
    model->searchClamours(pid);
}


void ClamourPage::on_getPetitionSupportButton_clicked()
{
    int nWindow = ui->petitionWindowSpinbox->value();
    model->getPetitionSupport(nWindow);
}

void ClamourPage::on_petitionSupportView_customContextMenuRequested(const QPoint &pos)
{
    QModelIndex index = ui->petitionSupportView->indexAt(pos);
    if (index.isValid())
        petitionViewContextMenu->exec(QCursor::pos());
}

void ClamourPage::searchHighlightedPetition()
{
    QModelIndexList indexes = ui->petitionSupportView->selectionModel()->selectedIndexes();
    if (indexes.size() >= 1)
    {
        QModelIndex pidIndex = indexes.at(0);
        QString pid = ui->petitionSupportView->model()->data(pidIndex).toString();
        ui->searchClamourEdit->setText(pid);
        ui->tabWidget->setCurrentIndex(1);
        ui->searchClamourButton->animateClick();
    }
}
